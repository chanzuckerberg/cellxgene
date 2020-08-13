import logging

import numpy as np
import pandas as pd


def get_dtype_of_array(array: pd.Series):
    return get_dtype_and_schema_of_array(array.dtype, array)[0]


def get_schema_type_hint_of_array(array: pd.Series):
    return get_dtype_and_schema_of_array(array.dtype, array)[1]


def get_dtype_and_schema_of_array(array: pd.Series):
    return (get_dtype_from_dtype(array.dtype, array), get_schema_type_hint_from_dtype(array.dtype))


def get_dtype_from_dtype(dtype, array_values=None):
    dtype_name = dtype.name
    dtype_kind = dtype.kind

    if dtype_name == np.float32 or dtype_name == np.int32:
        return dtype_name
    if dtype_name == "bool":
        return np.uint8
    if dtype_name == "object" and dtype_kind == "O":
        return np.unicode
    if dtype_name == "category":
        return get_dtype_from_dtype(dtype.categories)

    if can_cast_to_float32(dtype):
        return np.float32
    if can_cast_to_int32(dtype, array_values):
        return np.int32

    raise TypeError(f"Annotations of type {dtype} are unsupported.")


def get_schema_type_hint_from_dtype(dtype):
    if dtype == np.float32:
        return {"type": "float32"}
    if dtype == np.int32:
        return {"type": "int32"}
    if dtype == "bool":
        return {"type": "boolean"}
    if dtype == "object" and dtype.kind == "O":
        return {"type": "string"}
    if dtype == "category":
        return {"type": "categorical", "categories": dtype.categories.tolist()}

    raise TypeError(f"Annotations of type {dtype} are unsupported.")


def can_cast_to_float32(dtype):
    if dtype.kind == "f":
        if not np.can_cast(dtype, np.float32):
            logging.warning(f"Type {dtype.name} will be converted to 32 bit float and may lose precision.")
        return True
    return False


def can_cast_to_int32(dtype, array_values=None):
    if dtype.kind in ["i", "u"]:
        if np.can_cast(dtype, np.int32):
            return True
        ii32 = np.iinfo(np.int32)
        if array_values and (array_values.min() >= ii32.min and array_values.max() <= ii32.max) or array_values.empty:
            return True
    return False
