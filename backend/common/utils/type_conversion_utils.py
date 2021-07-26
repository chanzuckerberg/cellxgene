from typing import Union, Tuple
import logging

import numpy as np
import pandas as pd

"""
These routines drive all type inference for the schema generation and the
FBS (REST OTA) encoding. They are also used for CXG generation.


H5AD Type                       REST              REST
(ndarray, Series, Index)        FBS encoding      schema type       ERROR/exceptions
----------------------------    --------------    ---------------   ----------------------
bool_/bool                      uint8             boolean
(u)int8, (u)int16, int32        int32             int32
uint32, (u)it64                 int32             int32             CHECKS value bounds
float16, float32, float64       float32           float32[0]

categorical[T is numeric[4]]:
    hasna = False               T                 categorical[1]
    hasna = True                float32           categorical[1]    CHECKS value bounds

categorical[T not numeric]      JSON/str          categorical[1,2]

(other object)                  JSON/str          string

(all other)                                                         Always an ERROR[3]


Notes:
[0] IEEE format, includes non-finite numbers (NaN, Inf, ...)
[1] with NO categories enumerated (client side does it to handle rounding)
[2] NA (undefined) categories are assigned a JSON null value
[3] Includes all other numpy types:  datetime, complex, etc.
[4] means float, int, uint (dtype.kind in ['i','u','f'])

"""


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


def get_schema_type_hint_from_dtype(dtype) -> dict:
    res = _get_type_info_from_dtype(dtype)
    if res is None:
        raise TypeError(f"Annotations of type {dtype} are unsupported.")
    else:
        return res[1]


def _get_type_info_from_dtype(dtype) -> Union[Tuple[np.dtype, dict], None]:
    """
    Best-effort to determine encoding type and schema hint from a dtype.
    If this is not possible, or the type is unsupported, return None.
    """
    if dtype.kind == "b":
        return (np.uint8, {"type": "boolean"})

    if dtype.kind == "U":
        return (np.dtype(str), {"type": "string"})

    if dtype.kind in ["i", "u"]:
        if np.can_cast(dtype, np.int32):
            return (np.int32, {"type": "int32"})

    if dtype.kind == "f":
        if np.can_cast(dtype, np.float32):
            return (np.float32, {"type": "float32"})

    return None


def _get_type_info(array: Union[np.ndarray, pd.Series]) -> Tuple[np.dtype, dict]:
    """
    Determine encoding type and schema hint from an array.  This allows more
    flexible casting than may be possible by using just the dtype, as it can
    account for category types and array values.
    """
    if (
        not isinstance(array, np.ndarray)
        and not isinstance(array, pd.Series)
        and not isinstance(array, pd.Index)
        and not hasattr(array, "dtype")
    ):
        raise TypeError("Unsupported data type.")

    dtype = array.dtype

    res = _get_type_info_from_dtype(dtype)
    if res is not None:
        return res

    if dtype.kind == "O":
        if dtype.name == "category":
            # Sometimes CategoricalDType can be encoded as int or float without further fuss.
            # Do not specify the categories in the schema - let the client-side figure it out
            # on its own.  Utilize Series.to_numpy() to do casting that handles categorical
            # NA/NaN (missing or undefined) categories.
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

    if dtype.kind == "f" and can_cast_to_float32(array):
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
