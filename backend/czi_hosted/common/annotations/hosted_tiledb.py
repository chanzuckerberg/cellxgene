import json
import os
import re
import time

import pandas as pd
import tiledb
from flask import current_app

from backend.czi_hosted.common.annotations.annotations import Annotations
from backend.common_utils.errors import AnnotationCategoryNameError
from backend.czi_hosted.common.utils.sanitization_utils import sanitize_values_in_list
from backend.common_utils.type_conversion_utils import get_dtypes_and_schemas_of_dataframe, get_dtype_of_array
from backend.czi_hosted.db.cellxgene_orm import Annotation


class AnnotationsHostedTileDB(Annotations):
    CXG_ANNO_COLLECTION = "cxg_anno_collection"

    def __init__(self, directory_path, db):
        super().__init__()
        self.db = db
        if directory_path[-1] == "/":
            self.directory_path = directory_path
        else:
            self.directory_path = directory_path + "/"

    def check_category_names(self, df):
        original_category_names = df.keys().to_list()
        sanitized_category_names = set(sanitize_values_in_list(original_category_names).values())
        unsanitary_original_category_names = set(original_category_names).difference(sanitized_category_names)
        if unsanitary_original_category_names:
            raise AnnotationCategoryNameError(
                f"{unsanitary_original_category_names} are not valid category names, please resubmit"
            )

    def get_user_name(self):
        return current_app.auth.get_user_name()

    def get_user_id(self):
        return current_app.auth.get_user_id()

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
        user_id = self.get_user_id()
        if user_id is None:
            return
        dataset_name = data_adaptor.get_location()
        dataset_id = self.db.get_or_create_dataset(dataset_name)

        annotation_object = self.db.query_for_most_recent(
            Annotation, [Annotation.user_id == user_id, Annotation.dataset_id == dataset_id]
        )
        if annotation_object:
            if annotation_object.tiledb_uri == "":
                # this mean the user has removed all the categories.
                return None
            try:
                df = tiledb.open(annotation_object.tiledb_uri)
            except tiledb.TileDBError:
                # don't crash if the annotations file is missing or can't be read.
                current_app.logger.warning(f"Cannot read annotation file: {annotation_object.tiledb_uri}")
                return None
            pandas_df = self.convert_to_pandas_df(df, annotation_object.schema_hints)
            return pandas_df
        else:
            return None

    def convert_to_pandas_df(self, tileDBArray, schema_hints):
        repr_meta = None
        index_dims = None
        schema_hints = json.loads(schema_hints)

        if "__pandas_attribute_repr" in tileDBArray.meta:
            # backwards compatibility... unsure if necessary at this point
            repr_meta = json.loads(tileDBArray.meta["__pandas_attribute_repr"])
        if "__pandas_index_dims" in tileDBArray.meta:
            index_dims = json.loads(tileDBArray.meta["__pandas_index_dims"])

        data = tileDBArray[:]
        indexes = list()

        for col_name, col_val in data.items():
            # If the column values are byte literals, decode them
            if isinstance(col_val[0], bytes):
                col_val = [value.decode("utf-8") for value in col_val]

            if schema_hints and col_name in schema_hints:
                type = schema_hints.get(col_name).get("type")
                if type and type == "categorical":
                    new_col = pd.Series(col_val, dtype="category")
                    data[col_name] = new_col
            elif repr_meta and col_name in repr_meta:
                new_col = pd.Series(col_val, dtype=repr_meta[col_name])
                data[col_name] = new_col
            elif index_dims and col_name in index_dims:
                new_col = pd.Series(col_val, dtype=index_dims[col_name])
                data[col_name] = new_col
                indexes.append(col_name)

        new_df = pd.DataFrame.from_dict(data)
        if len(indexes) > 0:
            new_df.set_index(indexes, inplace=True)

        return new_df

    def write_labels(self, df, data_adaptor):
        auth_user_id = self.get_user_id()
        user_name = self.get_user_name()
        timestamp = time.time()
        dataset_location = data_adaptor.get_location()
        dataset_id = self.db.get_or_create_dataset(dataset_location)
        dataset_name = data_adaptor.get_title()
        user_id = self.db.get_or_create_user(auth_user_id)
        """
        NOTE: The uri contains the dataset name, user name and a timestamp as a convenience for debugging purposes.
        People may have the same name and time.time() can be server dependent.
        See - https://docs.python.org/2/library/time.html#time.time

        The annotations objects in the database should be used as the source of truth about who an annotation belongs
        to (for authorization purposes) and what time it was created (for garbage collection).
        """
        uri = f"{self.directory_path}{dataset_name}/{user_name}/{timestamp}"
        if uri.startswith("s3://"):
            pass
        else:
            os.makedirs(uri, exist_ok=True)
        _, dataframe_schema_type_hints = get_dtypes_and_schemas_of_dataframe(df)
        if not df.empty:
            self.check_category_names(df)
            # convert to tiledb datatypes

            for col in df:
                df[col] = df[col].astype(get_dtype_of_array(df[col]))
            tiledb.from_pandas(uri, df, sparse=True)
        else:
            uri = ""

        annotation = Annotation(
            tiledb_uri=uri,
            user_id=user_id,
            dataset_id=str(dataset_id),
            schema_hints=json.dumps(dataframe_schema_type_hints),
        )
        self.db.session.add(annotation)
        self.db.session.commit()

    def update_parameters(self, parameters, data_adaptor):
        params = {}
        params["annotations"] = True
        params["user_annotation_collection_name_enabled"] = False

        if self.ontology_data:
            params["annotations_cell_ontology_enabled"] = True
            params["annotations_cell_ontology_terms"] = self.ontology_data
        else:
            params["annotations_cell_ontology_enabled"] = False

        parameters.update(params)
