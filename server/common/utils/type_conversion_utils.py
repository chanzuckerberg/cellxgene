import logging

import numpy as np
import pandas as pd


def get_dtype_of_array(array: pd.Series):
    return get_dtype_and_schema_of_array(array)[0]


def get_schema_type_hint_of_array(array: pd.Series):
    return get_dtype_and_schema_of_array(array)[1]


def get_dtype_and_schema_of_array(array: pd.Series):
    return (get_dtype_from_dtype(array.dtype, array_values=array),
            get_schema_type_hint_from_dtype(array.dtype, array_values=array))


def get_dtype_from_dtype(dtype, array_values=None):
    """
    Given a data type, finds the equivalent data type that the array should be encoded as. Notably, this is relevant
    for 64 bit values which will get downcast to 32 bit.
    """

    dtype_name = dtype.name
    dtype_kind = dtype.kind

    if dtype == np.float32 or dtype == np.int32:
        return dtype
    if dtype_name == "bool":
        return np.uint8
    if dtype_name == "object" and dtype_kind == "O":
        return np.unicode
    if dtype_name == "category":
        return get_dtype_from_dtype(dtype.categories.dtype, dtype.categories)

    if can_cast_to_float32(dtype):
        return np.float32
    if can_cast_to_int32(dtype, array_values):
        return np.int32

    raise TypeError(f"Annotations of type {dtype} are unsupported.")


def get_schema_type_hint_from_dtype(dtype, array_values=None):
    """
    Returns a dictionary that contains type hints about the data type given, especially if the data type is 64 bit
    and will be downcast to 32 bit.
    """

    dtype_name = dtype.name
    dtype_kind = dtype.kind

    if dtype == np.float32 or dtype == np.int32:
        return {"type": dtype_name}
    if dtype_name == "bool":
        return {"type": "boolean"}
    if dtype_name == "object" and dtype_kind == "O":
        return {"type": "string"}
    if dtype_name == "category":
        return {"type": "categorical", "categories": dtype.categories.tolist()}

    if can_cast_to_float32(dtype):
        return {"type": "float32"}
    if can_cast_to_int32(dtype, array_values):
        return {"type": "int32"}

    raise TypeError(f"Annotations of type {dtype} are unsupported.")


def can_cast_to_float32(dtype):
    if dtype.kind == "f":
        if not np.can_cast(dtype, np.float32):
            logging.warning(f"Type {dtype.name} will be converted to 32 bit float and may lose precision.")
        return True
    return False


def can_cast_to_int32(dtype, array_values=None):
    """
    A type can be cast to 32 bit, overriding the numpy `cast_cast` function if the values in the array that are of
    the higher precision type has values that are entirely within the range of the downcast type.
    """

    if dtype.kind in ["i", "u"]:
        if np.can_cast(dtype, np.int32):
            return True
        ii32 = np.iinfo(np.int32)
        if not array_values.empty and (
                array_values.min() >= ii32.min and array_values.max() <= ii32.max) or array_values.empty:
            return True
    return False
