import logging

import numpy as np
import pandas as pd


def get_dtype_of_array(array: pd.Series):
    return get_dtype_and_schema_of_array(array)[0]


def get_schema_type_hint_of_array(array: pd.Series):
    return get_dtype_and_schema_of_array(array)[1]


def get_dtype_and_schema_of_array(array: pd.Series):
    dtype = array.dtype
    data_kind = dtype.kind

    if dtype == np.float32:
        return (dtype, {"type": "float32"})
    if dtype == np.int32:
        return (dtype, {"type": "int32"})
    if dtype == np.bool_:
        return (np.uint8, {"type": "boolean"})
    if dtype == np.str:
        return (np.unicode, {"type": "string"})
    if dtype == "category":
        type, _ = get_dtype_and_schema_of_array(dtype.categories)
        return (type, {"type": "categorical", "categories": dtype.categories.tolist()})

    # Check whether the data type can be cast to either float32 or int32
    if can_cast_to_float32(array):
        return (np.float32, {"type": "float32"})
    if can_cast_to_int32(array):
        return (np.int32, {"type": "int32"})

    # Check whether the data is an object
    if data_kind == "O" and dtype == "object":
        return (np.unicode, {"type": "string"})

    raise TypeError(f"Annotations of type {dtype} are unsupported.")


def can_cast_to_float32(array: pd.Series):
    if array.dtype.kind == "f":
        if not np.can_cast(array.dtype, np.float32):
            logging.warning(f"Array {array.name} will be converted to 32 bit float and may lose precision.")
        return True
    return False


def can_cast_to_int32(array: pd.Series):
    if array.dtype.kind in ["i", "u"]:
        if np.can_cast(array.dtype, np.int32):
            return True
        ii32 = np.iinfo(np.int32)
        if array.min() >= ii32.min and array.max() <= ii32.max:
            return True
    return False
