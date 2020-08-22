import json
import os
import re
import time

import numpy as np
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
            pandas_df = self.convert_to_pandas_df(df, annotation_object.schema_hints)
            return pandas_df
        else:
            return None

    def convert_to_pandas_df(self, tileDBArray, json_schema_hints):
        values = tileDBArray[:]
        schema_hints = json.loads(json_schema_hints)
        index = None

        dataframe_data = {}

        for column_name, column_values in values.items():
            print(column_name)
            if column_name == 'index':
                index = pd.Series(column_values, dtype=np.unicode, name=column_name)
                index = index.astype(np.unicode)
            else:
                type_hint = schema_hints.get(column_name)
                if type_hint:
                    type = type_hint.get("type")
                    if type == "boolean":
                        series = pd.Series(column_values, dtype=np.bool_)
                    elif type == "float32":
                        series = pd.Series(column_values, dtype=np.float32)
                    elif type == "int32":
                        series = pd.Series(column_values, dtype=np.int32)
                    elif type == "categorical":
                        series = pd.Series(column_values, dtype='category')
                    else:
                        series = pd.Series(column_values, dtype=type)
                else:
                    series = pd.Series(column_values, dtype=column_values.dtype, name=column_name)

                dataframe_data[column_name] = series

        dataframe = pd.DataFrame(data=dataframe_data)
        if index is not None:
            dataframe = dataframe.set_index(index)
        return dataframe

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
