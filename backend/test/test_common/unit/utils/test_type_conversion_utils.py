import unittest
import logging
from parameterized import parameterized_class

import numpy as np
import pandas as pd
from pandas import Series, DataFrame
from scipy import sparse

from backend.common.utils.type_conversion_utils import (
    get_encoding_dtype_of_array,
    get_schema_type_hint_of_array,
    get_dtypes_and_schemas_of_dataframe,
    get_dtype_and_schema_of_array,
    get_schema_type_hint_from_dtype,
)


class TestTypeConversionUtils(unittest.TestCase):
    def test__get_dtypes_and_schemas_of_dataframe__dtype_and_schema_returns_as_expected(self):
        float_array = Series(data=[1, 2, 3], dtype=np.dtype(np.float64))
        category_array = Series(data=["a", "b", "b"], dtype="category")
        dataframe = DataFrame({"float_array": float_array, "category_array": category_array})

        expected_data_types_dict = {"float_array": np.float32, "category_array": str}
        expected_schema_type_hints_dict = {
            "float_array": {"type": "float32"},
            "category_array": {"type": "categorical", "categories": ["a", "b"]},
        }

        actual_dataframe_data_types, actual_dataframe_schema_type_hints = get_dtypes_and_schemas_of_dataframe(dataframe)

        self.assertEqual(expected_data_types_dict, actual_dataframe_data_types)
        self.assertEqual(expected_schema_type_hints_dict, actual_dataframe_schema_type_hints)

    def test__get_schema_type_hint_from_dtype(self):
        self.assertEqual(get_schema_type_hint_from_dtype(np.dtype(np.bool_)), {"type": "boolean"})

        for dtype in [np.int8, np.int8, np.int16, np.uint16, np.int32]:
            self.assertEqual(get_schema_type_hint_from_dtype(np.dtype(dtype)), {"type": "int32"})
        for dtype in [np.uint32, np.int64, np.uint64]:
            with self.assertRaises(TypeError):
                get_schema_type_hint_from_dtype(np.dtype(dtype))

        for dtype in [np.float16, np.float32, np.float64]:
            self.assertEqual(get_schema_type_hint_from_dtype(np.dtype(dtype)), {"type": "float32"})

        for dtype in [np.dtype(object), np.dtype(str)]:
            self.assertEqual(get_schema_type_hint_from_dtype(dtype), {"type": "string"})


# Credit: https://stackoverflow.com/questions/35871815/python-3-unit-testing-assert-logger-not-called/64774103#64774103
class AssertNoLog:
    def assertNoLogs(self, logger, level):
        """functions as a context manager.  To be introduced in python 3.10"""

        class AssertNoLogsContext(unittest.TestCase):
            def __init__(self, logger, level):
                self.logger = logger
                self.level = level
                self.context = self.assertLogs(logger, level)

            def __enter__(self):
                """enter self.assertLogs as context manager, and log something"""
                self.initial_logmsg = "sole message"
                self.cm = self.context.__enter__()
                self.logger.log(self.level, self.initial_logmsg)
                return self.cm

            def __exit__(self, exc_type, exc_val, exc_tb):
                """cleanup logs, and then check nothing extra was logged"""
                # assertLogs.__exit__ should never fail because of initial msg
                self.context.__exit__(exc_type, exc_val, exc_tb)
                if len(self.cm.output) > 1:
                    """override any exception passed to __exit__"""
                    self.context._raiseFailure(
                        "logs of level {} or higher triggered on {} : {}".format(
                            logging.getLevelName(self.level), self.logger.name, self.cm.output[1:]
                        )
                    )

        return AssertNoLogsContext(logger, level)


"""
See table of expected cases in type_conversion_utils.py.

This probes all edge cases. Each case is a dict containing keys:
    - data - the array to be introspected
    - throws - if not None, the expected Error (eg, TypeError)
    - expected_encoding_dtype - upon success
    - expected_schema_hint - upon success
    - logs - if not None, specify expected log output
"""

bool_OK_cases = [
    {
        "data": data,
        "expected_encoding_dtype": np.uint8,
        "expected_schema_hint": {"type": "boolean"},
    }
    for data in [
        np.array([0, 1, 0, 1], dtype=np.bool_),
        pd.Series(np.array([0, 1, 0, 1], dtype=np.bool_)),
        # pd.Index with bools doesn't really make any sense...and becomes dtype=object
    ]
]

int_OK_cases = [
    {
        "data": data,
        "expected_encoding_dtype": np.int32,
        "expected_schema_hint": {"type": "int32"},
    }
    for dtype in [np.int8, np.uint8, np.int16, np.uint16, np.int32, np.uint32, np.int64, np.uint64]
    for data in [
        np.arange(0, 1000, dtype=dtype),
        pd.Series(np.arange(0, 1000, dtype=dtype)),
        pd.Index(np.arange(0, 1000, dtype=dtype)),
        sparse.csr_matrix((10, 100), dtype=dtype),
    ]
]

float_OK_cases = [
    {
        "data": data,
        "expected_encoding_dtype": np.float32,
        "expected_schema_hint": {"type": "float32"},
        "logs": None if data.dtype != np.float64 else {"level": logging.WARNING, "output": "may lose precision"},
    }
    for dtype in [np.float16, np.float32, np.float64]
    for data in [
        np.arange(-128, 1000, dtype=dtype),
        pd.Series(np.arange(-128, 1000, dtype=dtype)),
        pd.Index(np.arange(-129, 1000, dtype=dtype)),
        np.array([-np.nan, np.NINF, -1, np.NZERO, 0, np.PZERO, 1, np.PINF, np.nan], dtype=dtype),
        np.array([np.finfo(dtype).min, 0, np.finfo(dtype).max], dtype=dtype),
        sparse.csr_matrix((10, 100), dtype=dtype),
    ]
]


numeric_ERR_cases = [
    {
        "data": data,
        "throws": TypeError,
    }
    for data in [
        np.array([np.iinfo(np.int64).min, np.iinfo(np.int64).max], dtype=np.int64),
        np.array([np.iinfo(np.uint64).min, np.iinfo(np.uint64).max], dtype=np.uint64),
        np.array([np.iinfo(np.uint32).min, np.iinfo(np.uint32).max], dtype=np.uint32),
    ]
]


string_OK_cases = [
    {
        "data": data,
        "expected_encoding_dtype": np.dtype(str),
        "expected_schema_hint": {"type": "string"},
    }
    for data in [
        np.array(["a", "b", "c"]),
        np.array(["a", "b", "c"], dtype="object"),
        pd.Series(["a", "b", "c"]),
        pd.Index(["a", "b", "c"]),
        np.array(["a", [], {}, None, True, False, 383.2], dtype="object"),
    ]
]

category_nonnumeric_OK_cases = [
    {
        "data": data,
        "expected_encoding_dtype": np.dtype(str),
        "expected_schema_hint": {"type": "categorical", "categories": data.dtype.categories.to_list()},
    }
    for data in [
        pd.Series(["a", "b", "c"], dtype="category"),
        pd.Series(["a", "b", "c", 0, 1, 2], dtype="category"),
        pd.Series(["a", "b", "c"], dtype="category").cat.remove_categories(["b"]),
        pd.Series(["a", "b", "c", 0, 1, 2], dtype="category").cat.remove_categories(["b", 0]),
    ]
]

category_numeric_OK_cases = [
    # numeric, no NA/NaN, int
    *[
        {
            "data": data,
            "expected_encoding_dtype": np.int32,
            "expected_schema_hint": {"type": "categorical"},
        }
        for dtype in [np.int8, np.uint8, np.int16, np.uint16, np.int32, np.uint32, np.int64, np.uint64]
        for data in [
            pd.Series(np.array([0, 1, 2], dtype=dtype), dtype="category"),
        ]
    ],
    # numeric, no NA/NaN, float
    *[
        {
            "data": data,
            "expected_encoding_dtype": np.float32,
            "expected_schema_hint": {"type": "categorical"},
            "logs": {"level": logging.WARNING, "output": "may lose precision"},
        }
        for dtype in [np.float16, np.float32, np.float64]
        for data in [
            pd.Series(np.array([0, 1, 2], dtype=dtype), dtype="category"),
            pd.Series(np.array([0, 1, 2], dtype=dtype), dtype="category").cat.remove_categories([1]),
            pd.Categorical(np.array([0, 1, 2], dtype=dtype)),
        ]
    ],
    # numeric, has NA-induced cast to float32
    *[
        {
            "data": data,
            "expected_encoding_dtype": np.float32,
            "expected_schema_hint": {"type": "categorical"},
            "logs": {"level": logging.WARNING, "output": "may lose precision"},
        }
        for dtype in [
            np.int8,
            np.uint8,
            np.int16,
            np.uint16,
            np.int32,
            np.uint32,
            np.int64,
            np.uint64,
            np.float16,
            np.float32,
            np.float64,
        ]
        for data in [
            pd.Series(np.array([0, 1, 2], dtype=dtype), dtype="category").cat.remove_categories([1]),
            pd.Categorical(np.array([0, 1, 2], dtype=dtype), categories=np.array([0, 1], dtype=dtype)),
        ]
    ],
]

category_ERR_cases = [
    # catch expected categorical exceptions for Int64(etc) that have large values
    {
        "data": data,
        "throws": TypeError,
    }
    for data in [
        pd.Categorical(np.array([np.iinfo(np.int64).min, np.iinfo(np.int64).max], dtype=np.int64)),
        pd.Categorical(np.array([np.iinfo(np.uint64).min, np.iinfo(np.uint64).max], dtype=np.uint64)),
        pd.Categorical(np.array([np.iinfo(np.uint32).min, np.iinfo(np.uint32).max], dtype=np.uint32)),
    ]
]

object_OK_cases = [
    {
        "data": data,
        "expected_encoding_dtype": np.dtype(str),
        "expected_schema_hint": {"type": "string"},
    }
    for data in [
        np.array(["a", True, 1, [], {}], dtype="object"),
        pd.Series(["a", True, 1, [], {}], dtype="object"),
        pd.Index(["a", True, 1, [], {}], dtype="object"),
    ]
]

err_cases = [
    {"data": np.array, "throws": TypeError}
    for data in [
        np.ones((10,), dtype=np.complex64),
        np.ones((10,), dtype=np.complex128),
        np.array([b"foobar"], dtype=np.bytes_),
        np.ones((10,), dtype=np.void),
        np.arange("2005-02", "2005-03", dtype="datetime64[D]"),
        np.arange("2005-02", "2005-03", dtype="datetime64[D]") - np.datetime64("2008-01-01"),
        [],
        {},
    ]
]

test_cases = [
    *bool_OK_cases,
    *int_OK_cases,
    *float_OK_cases,
    *numeric_ERR_cases,
    *string_OK_cases,
    *category_nonnumeric_OK_cases,
    *category_numeric_OK_cases,
    *category_ERR_cases,
    *object_OK_cases,
    *err_cases,
]


@parameterized_class(test_cases)
class TestTypeInference(unittest.TestCase, AssertNoLog):
    def test_type_inference(self):
        throws = getattr(self, "throws", None)
        if throws:
            with self.assertRaises(throws):
                get_dtype_and_schema_of_array(self.data)
            with self.assertRaises(throws):
                get_encoding_dtype_of_array(self.data)
            with self.assertRaises(throws):
                get_schema_type_hint_of_array(self.data)

        else:
            logs = getattr(self, "logs", None)
            if logs is not None:
                with self.assertLogs(level=logs["level"]) as logger:
                    encoding_dtype, schema_hint = get_dtype_and_schema_of_array(self.data)
                    self.assertEqual(encoding_dtype, self.expected_encoding_dtype)
                    self.assertEqual(schema_hint, self.expected_schema_hint)
                    self.assertIn(logs["output"], logger.output[0])

            else:
                with self.assertNoLogs(logging.getLogger(), logging.WARNING):
                    encoding_dtype, schema_hint = get_dtype_and_schema_of_array(self.data)
                    self.assertEqual(encoding_dtype, self.expected_encoding_dtype)
                    self.assertEqual(schema_hint, self.expected_schema_hint)

            # also test the other public API
            self.assertEqual(get_encoding_dtype_of_array(self.data), self.expected_encoding_dtype)
            self.assertEqual(get_schema_type_hint_of_array(self.data), self.expected_schema_hint)
