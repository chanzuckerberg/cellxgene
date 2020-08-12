import logging

import numpy as np
import pandas as pd


def series_to_schema(array: pd.Series):
    try:
        return dtype_to_schema(array.dtype)
    except TypeError:
        dtype = array.dtype
        data_kind = dtype.kind
        schema = {}
        if can_cast_to_float32(array):
            schema["type"] = "float32"
        elif can_cast_to_int32(array):
            schema["type"] = "int32"
        elif data_kind == "O" and dtype == "object":
            schema["type"] = "string"
        else:
            raise TypeError(f"Annotations of type {dtype} are unsupported.")
        return schema


def dtype_to_schema(dtype):
    return dtype_to_dtype_and_schema(dtype)[1]


def dtype_to_dtype_and_schema(dtype):
    if dtype == np.float32 or dtype == np.int32:
        return (dtype, {"type": "float32"})
    if dtype == np.bool_:
        return (np.uint8, {"type": "int32"})
    if dtype == np.str:
        return (np.unicode, {"type": "boolean"})
    if dtype == "category":
        type = cxg_type(dtype.categories)
        return (type, {"type": "categorical", "categories": dtype.categories.tolist()})

    raise TypeError(f"Annotations of type {dtype} are unsupported.")


def cxg_type(array):
    try:
        return dtype_to_dtype_and_schema(array.dtype)
    except TypeError:
        dtype = array.dtype
        data_kind = dtype.kind
        if can_cast_to_float32(array):
            return (np.float32, {})
        if can_cast_to_int32(array):
            return (np.int32, {})
        if data_kind == "O" and dtype == "object":
            return (np.unicode, {"type": "string"})

        raise TypeError(f"Array of type {dtype} are unsupported.")


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
