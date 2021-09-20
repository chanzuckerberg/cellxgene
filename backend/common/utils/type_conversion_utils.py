from typing import Union, Tuple
import logging

import numpy as np
import pandas as pd

"""
These routines drive all type inference for the schema generation and the
FBS (REST OTA) encoding.


H5AD Type                       REST              REST
(ndarray, Series, Index)        FBS encoding      schema type       ERROR/exceptions
----------------------------    --------------    ---------------   ----------------------
bool_/bool                      uint8             boolean
(u)int8, (u)int16, int32        int32             int32
uint32, (u)int64                int32             int32             CHECKS value bounds
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


def get_encoding_dtype_of_array(array: Union[np.ndarray, pd.Series, pd.Index]) -> np.dtype:
    return _get_type_info(array)[0]


def get_schema_type_hint_of_array(array: Union[np.ndarray, pd.Series, pd.Index]) -> dict:
    return _get_type_info(array)[1]


def get_dtype_and_schema_of_array(array: Union[np.ndarray, pd.Series, pd.Index]) -> Tuple[np.dtype, dict]:
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

    This should be a subset of the cases which are supported by
    _get_type_info().  The latter should be preferred if the array (values)
    are available for typing.
    """
    if dtype.kind == "b":
        return (np.uint8, {"type": "boolean"})

    if dtype.kind == "U":
        return (np.dtype(str), {"type": "string"})

    if dtype.kind in ["i", "u"]:
        if np.can_cast(dtype, np.int32):
            return (np.int32, {"type": "int32"})

    if dtype.kind == "f":
        _float64_warning(dtype)
        return (np.float32, {"type": "float32"})

    if dtype.kind == "O" and not dtype.name == "category":
        return (np.dtype(str), {"type": "string"})

    return None


def _get_type_info(array: Union[np.ndarray, pd.Series, pd.Index]) -> Tuple[np.dtype, dict]:
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

    if dtype.kind in ["i", "u"] and _can_cast_array_values_to_int32(array):
        return (np.int32, {"type": "int32"})

    if dtype.kind == "f":
        _float64_warning(array.dtype)
        return (np.float32, {"type": "float32"})

    raise TypeError(f"Annotations of type {dtype} are unsupported.")


def _float64_warning(dtype):
    """
    Warn the user if we are down-casting a float64 to float32, and may potentially lose information.
    """
    if dtype.kind == "f" and not np.can_cast(dtype, np.float32):
        logging.warning(f"Type {dtype.name} will be converted to 32 bit float and may lose precision.")


def _can_cast_array_values_to_int32(array: Union[np.ndarray, pd.Series, pd.Index]) -> bool:
    """
    Return true if the (U)INT array values can be safely cast to int32.  We allow size reducing
    casts (ie, int64 to int32) if no actual values require the larger size (ie, actual values
    can be represented by the smaller type).
    """
    assert array.dtype.kind in ["u", "i"]

    if np.can_cast(array.dtype, np.int32):
        return True

    if array.size == 0:
        return True

    int32_machine_limits = np.iinfo(np.int32)
    if array.min() >= int32_machine_limits.min and array.max() <= int32_machine_limits.max:
        return True

    return False


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
