from typing import Union, Tuple
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


def get_encoding_dtype_of_array(array: Union[np.ndarray, pd.Series]) -> np.dtype:
    return _get_type_info(array)[0]


def get_schema_type_hint_of_array(array: Union[np.ndarray, pd.Series]) -> dict:
    return _get_type_info(array)[1]


def get_dtype_and_schema_of_array(array: Union[np.ndarray, pd.Series]) -> Tuple[np.dtype, dict]:
    """Return tuple (encoding_dtype, schema_type_hint)"""
    return _get_type_info(array)


def _get_type_info(array: Union[np.ndarray, pd.Series]) -> Tuple[np.dtype, dict]:
    """
    Given a data type, return the type information used to:
        * FBS encoding dtype: encode Matrix FBS format using this primitive type.
        * Client schema info: which will be passed to the client
    Only return information for those types that we have full support for in
    both the backend AND front-end.

    Note: in the front-end:
    * the front-end only supports int32 and float32 numerics

    In our FBS encoding:
    * pandas.Series are cast to arrays using Series.to_numpy(), with DEFAULT
      pandas type conversion rules (this handles, for example, categoricals
      with missing values)
    * only unicode strings are supported (no bytes or bytearray)
    * no arbitrary object support

    CategoricalDType are special, in that the schema hint wants to know about the
    category and categories, and we want to encode numeric categoricals
    """
    dtype = array.dtype

    if dtype.kind == "b" or dtype.name == "bool":
        return (np.uint8, {"type": "boolean"})

    if dtype.kind == "U":
        return (np.dtype(str), {"type": "string"})

    if dtype.kind == "O":
        if dtype.name == "category":
            assert isinstance(array, pd.Series)
            # Sometimes CategoricalDType can be encoded as int or float without further fuss.
            # Do not specify the categories in the schema - let the client-side figure it out
            # on its own.  Utilize Series.to_numpy() to
            if dtype.categories.dtype.kind in ["f", "i", "u"]:
                return (
                    _get_type_info(array.to_numpy())[0],
                    {"type": "categorical"},
                )
            else:
                return (np.dtype(str), {"type": "categorical", "categories": dtype.categories.to_list()})

        # all other extension types are str-encoded
        return (np.dtype(str), {"type": "string"})

    if dtype.kind in ["i", "u"] and can_cast_to_int(array, np.int32):
        return (np.int32, {"type": "int32"})

    if can_cast_to_float32(array):
        return (np.float32, {"type": "float32"})

    raise TypeError(f"Annotations of type {dtype} are unsupported.")


def can_cast_to_float32(array: Union[np.ndarray, pd.Series]) -> bool:
    """
    Returns True if the underlying type is float, of any width.
    Else, returns False
    """
    # Accept any float as long as magnitude can still be represented.
    if array.dtype.kind == "f":
        if not np.can_cast(array.dtype, np.float32):
            logging.warning(f"Type {array.dtype.name} will be converted to 32 bit float and may lose precision.")
        return True

    return False


def can_cast_to_int(array: Union[np.ndarray, pd.Series], target: np.dtype) -> bool:
    """
    Return true if the dtype can be safely cast to int dtype.  We allow size reducing casts
    (ie, int64 to int32) if no actual values require the larger size (ie, actual values
    can be represented by the smaller type).
    """
    if array.dtype.kind in ["u", "i"]:
        if np.can_cast(array.dtype, target):
            return True

        if array.size == 0:
            return True

        ii = np.iinfo(target)
        if array.min() >= ii.min and array.max() <= ii.max:
            return True

    return False


def convert_pandas_series_to_numpy(series_to_convert: pd.Series, dtype):
    if series_to_convert.hasnans and dtype == np.int32:
        logging.error("Cannot convert a pandas Series object to an integer dtype if it contains NA values.")

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
