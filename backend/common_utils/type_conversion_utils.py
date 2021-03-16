import logging

import numpy as np
import pandas as pd


def get_dtypes_and_schemas_of_dataframe(dataframe: pd.DataFrame):
    dtypes_by_column_name = {}
    schema_type_hints_by_column_name = {}

    for column_name, column_values in dataframe.items():
        (
            dtypes_by_column_name[column_name],
            schema_type_hints_by_column_name[column_name],
        ) = get_dtype_and_schema_of_array(column_values)

    return dtypes_by_column_name, schema_type_hints_by_column_name


def get_dtype_of_array(array: pd.Series):
    return get_dtype_and_schema_of_array(array)[0]


def get_schema_type_hint_of_array(array: pd.Series):
    return get_dtype_and_schema_of_array(array)[1]


def get_dtype_and_schema_of_array(array: pd.Series):
    return (
        get_dtype_from_dtype(array.dtype, array_values=array),
        get_schema_type_hint_from_dtype(array.dtype, array_values=array),
    )


def get_dtype_from_dtype(dtype, array_values=None):
    """
    Given a data type, finds the equivalent data type that the array should be encoded as. Notably, this is relevant
    for 64 bit values which will get downcast to 32 bit.
    """

    dtype_name = dtype.name
    dtype_kind = dtype.kind

    if dtype_name == "bool":
        return np.uint8
    if dtype_name == "object" and dtype_kind == "O":
        return str
    if dtype_name == "category":
        return get_dtype_from_dtype(dtype.categories.dtype, array_values)

    if can_cast_to_int32(dtype, array_values):
        return np.int32
    if can_cast_to_float32(dtype, array_values):
        return np.float32
    if not can_cast_to_float32(dtype, array_values):
        return np.float64

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

    if can_cast_to_int32(dtype, array_values):
        return {"type": "int32"}
    if can_cast_to_float32(dtype, array_values):
        return {"type": "float32"}
    if dtype_kind == "f" and not can_cast_to_float32(dtype, array_values):
        return {"type": "float64"}

    raise TypeError(f"Annotations of type {dtype} are unsupported.")


def can_cast_to_float32(dtype, array_values):
    """
    Optimistically returns True signifying that a type downcast to float32 is possible whenever the incoming type is
    a float.

    We also handle a special case here where the array is a Series object with integer categorical values AND NaNs.
    Since NaNs are floating points in numpy, we upcast the integer array to float32 and return True.
    """

    if dtype.kind == "f":
        if not np.can_cast(dtype, np.float32):
            logging.warning(f"Type {dtype.name} will be converted to 32 bit float and may lose precision.")

        return True

    if dtype.kind == "O" and array_values.hasnans:
        return True

    return False


def can_cast_to_int32(dtype, array_values=None):
    """
    A type can be cast to 32 bit, overriding the numpy `cast_cast` function if the values in the array that are of
    the higher precision type has values that are entirely within the range of the downcast type.
    """

    # Since a NaN is technically a float, any array that contains NaNs cannot be cast to an integer so immediately
    # return False.
    if array_values.hasnans:
        return False

    # If the array is categorical, then we need to order the array values so that functions min and max that occur
    # later, can function. They do not function on unordered categories.
    ordered_array_values = array_values
    if array_values.dtype.name == "category" and not array_values.cat.ordered:
        ordered_array_values = array_values.cat.as_ordered()

    if dtype.kind in ["i", "u"]:
        if np.can_cast(dtype, np.int32):
            return True
        ii32 = np.iinfo(np.int32)
        if (
            not ordered_array_values.empty
            and (ordered_array_values.min() >= ii32.min and ordered_array_values.max() <= ii32.max)
            or ordered_array_values.empty
        ):
            return True
    return False


def convert_pandas_series_to_numpy(series_to_convert: pd.Series, dtype):
    if series_to_convert.hasnans and dtype == np.int32:
        logging.error("Cannot convert a pandas Series object to an integer dtype if it contains NaNs.")

    return series_to_convert.to_numpy(dtype)


def convert_string_to_value(value: str):
    """convert a string to value with the most appropriate type"""
    if value.lower() == "true":
        return True
    if value.lower() == "false":
        return False
    if value == "null":
        return None
    try:
        return eval(value)
    except:  # noqa E722
        return value
