import json
import os
import re
import time

import pandas as pd
import tiledb
from flask import current_app

from server.common.annotations.annotations import Annotations
from server.common.errors import AnnotationCategoryNameError
from server.common.utils.sanitization_utils import sanitize_values_in_list
from server.common.utils.type_conversion_utils import get_dtypes_and_schemas_of_dataframe, get_dtype_of_array
from server.db.cellxgene_orm import CellxGeneDataset, Annotation


class AnnotationsHostedTileDB(Annotations):
    CXG_ANNO_COLLECTION = "cxg_anno_collection"

    def __init__(self, directory_path, db):
        super().__init__()
        self.db = db
        self.directory_path = directory_path

    def check_category_names(self, df):
        original_category_names = df.keys().to_list()
        sanitized_category_names = set(sanitize_values_in_list(original_category_names).values())
        unsanitary_original_category_names = set(original_category_names).difference(sanitized_category_names)
        if unsanitary_original_category_names:
            raise AnnotationCategoryNameError(
                f"{unsanitary_original_category_names} are not valid category names, please resubmit")

    def is_safe_collection_name(self, name):
        """
        return true if this is a safe collection name
        this is ultra conservative. If we want to allow full legal file name syntax,
        we could look at modules like `pathvalidate`
        """
        if name is None:
            return False
        return re.match(r"^[\w\-]+$", name) is not None

    def set_collection(self, name):
        self.CXG_ANNO_COLLECTION = name

    def read_labels(self, data_adaptor):
        user_id = current_app.auth.get_user_id()
        dataset_name = data_adaptor.get_location()
        dataset_id = str(self.db.query(
            table_args=[CellxGeneDataset],
            filter_args=[CellxGeneDataset.name == dataset_name]
        )[0].id)

        annotation_object = self.db.query_for_most_recent(
            Annotation, [Annotation.user_id == user_id, Annotation.dataset_id == dataset_id]
        )
        if annotation_object:
            df = tiledb.open(annotation_object.tiledb_uri)
            pandas_df = self.convert_to_pandas_df(df)
            return pandas_df
        else:
            return None

    def convert_to_pandas_df(self, tileDBArray):
        repr_meta = None
        index_dims = None
        if '__pandas_attribute_repr' in tileDBArray.meta:
            # backwards compatibility... unsure if necessary at this point
            repr_meta = json.loads(tileDBArray.meta['__pandas_attribute_repr'])
        if '__pandas_index_dims' in tileDBArray.meta:
            index_dims = json.loads(tileDBArray.meta['__pandas_index_dims'])

        data = tileDBArray[:]
        indexes = list()

        for col_name, col_val in data.items():
            if repr_meta and col_name in repr_meta:
                new_col = pd.Series(col_val, dtype=repr_meta[col_name])
                data[col_name] = new_col
            elif index_dims and col_name in index_dims:
                new_col = pd.Series(col_val)
                data[col_name] = new_col
                indexes.append(col_name)

        new_df = pd.DataFrame.from_dict(data)
        if len(indexes) > 0:
            new_df.set_index(indexes, inplace=True)

        return new_df

    def write_labels(self, df, data_adaptor):

        user_id = current_app.auth.get_user_id()
        timestamp = time.time()
        dataset_name = data_adaptor.get_location()
        dataset_id = self.db.get_or_create_dataset(dataset_name)
        user_id = self.db.get_or_create_user(user_id)

        uri = f"{self.directory_path}-{dataset_name}-{user_id}-{timestamp}"
        if uri.startswith("s3://"):
            pass
        else:
            os.makedirs(uri, exist_ok=True)
        _, dataframe_schema_type_hints = get_dtypes_and_schemas_of_dataframe(df)

        annotation = Annotation(
            tiledb_uri=uri,
            user_id=user_id,
            dataset_id=str(dataset_id),
            schema_hints=json.dumps(dataframe_schema_type_hints)
        )
        if not df.empty:
            self.check_category_names(df)
            # convert to tiledb datatypes
            for col in df:
                df[col] = df[col].astype(get_dtype_of_array(df[col]))
            tiledb.from_pandas(uri, df)

        self.db.session.add(annotation)
        self.db.session.commit()
