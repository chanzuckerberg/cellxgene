import unittest
from time import time
from unittest.mock import patch

import numpy as np
from pandas import Series, DataFrame

from server.common.utils.type_conversion_utils import (
    can_cast_to_float32,
    can_cast_to_int32,
    get_dtype_of_array,
    get_schema_type_hint_of_array,
    get_dtypes_and_schemas_of_dataframe,
    convert_pandas_series_to_numpy,
)


class TestTypeConversionUtils(unittest.TestCase):
    def test__can_cast_to_float32__string_is_false(self):
        array_to_convert = Series(data=["1", "2", "3"], dtype=str)

        can_cast = can_cast_to_float32(array_to_convert.dtype, array_to_convert)

        self.assertFalse(can_cast)

    def test__can_cast_to_float32__float64_is_true_warning_outputted(self):
        array_to_convert = Series(data=[1, 2, 3], dtype=np.dtype(np.float64))

        with self.assertLogs(level="WARN") as logger:
            can_cast = can_cast_to_float32(array_to_convert.dtype, array_to_convert)
            self.assertIn("may lose precision", logger.output[0])

        self.assertTrue(can_cast)

    @patch("logging.warning")
    def test__can_cast_to_float32__float32_is_false(self, mock_log_warning):
        array_to_convert = Series(data=[1, 2, 3], dtype=np.dtype(np.float32))

        can_cast = can_cast_to_float32(array_to_convert.dtype, array_to_convert)

        self.assertTrue(can_cast)
        assert not mock_log_warning.called

    def test__can_cast_to_float32__categorical_float64_is_false(self):
        array_to_convert = Series(data=[1.1, 2.2, 3.3], dtype="category")

        can_cast = can_cast_to_float32(array_to_convert.dtype, array_to_convert)

        self.assertFalse(can_cast)

    def test__can_cast_to_float32__categorical_int64_with_nans_is_true(self):
        array_to_convert = Series(data=[1, 2, np.NaN], dtype="category")

        can_cast = can_cast_to_float32(array_to_convert.dtype, array_to_convert)

        self.assertTrue(can_cast)

    def test__can_cast_to_float_32__float_32_with_nans_is_true(self):
        array_to_convert = Series(data=[1, 2, np.NaN], dtype=np.dtype(np.float32))

        can_cast = can_cast_to_float32(array_to_convert.dtype, array_to_convert)

        self.assertTrue(can_cast)

    def test__can_cast_to_int32__string_is_false(self):
        array_to_convert = Series(data=["1", "2", "3"], dtype=str)

        can_cast = can_cast_to_int32(array_to_convert.dtype, array_to_convert)

        self.assertFalse(can_cast)

    def test__can_cast_to_int32__int64_is_true(self):
        array_to_convert = Series(data=["1", "2", "3"], dtype=np.dtype(np.int64))

        can_cast = can_cast_to_int32(array_to_convert.dtype, array_to_convert)

        self.assertTrue(can_cast)

    def test__can_cast_to_int32__int16_is_true(self):
        array_to_convert = Series(data=["1", "2", "3"], dtype=np.dtype(np.int16))

        can_cast = can_cast_to_int32(array_to_convert.dtype, array_to_convert)

        self.assertTrue(can_cast)

    def test__can_cast_to_int32__int64_with_large_value_is_false(self):
        array_to_convert = Series(data=["3000000000", "2", "3"], dtype=np.dtype(np.int64))

        can_cast = can_cast_to_int32(array_to_convert.dtype, array_to_convert)

        self.assertFalse(can_cast)

    def test__can_cast_to_int32__int64_with_nans_is_false(self):
        array_to_convert = Series(data=[np.NaN, "2", "3"], dtype="category")

        can_cast = can_cast_to_int32(array_to_convert.dtype, array_to_convert)

        self.assertFalse(can_cast)

    def test__get_dtype_of_array__supported_dtypes_return_as_expected(self):
        types = [np.float32, np.int32, np.bool_, str]
        expected_dtypes = [np.float32, np.int32, np.uint8, np.unicode]

        for test_type_index in range(len(types)):
            with self.subTest(
                f"Testing get_dtype_of_array with type {types[test_type_index].__name__}", i=test_type_index
            ):
                array = Series(data=[], dtype=types[test_type_index])
                self.assertEqual(get_dtype_of_array(array), expected_dtypes[test_type_index])

    def test__get_dtype_of_array__categories_return_as_expected(self):
        array = Series(data=["a", "b", "c"], dtype="category")
        expected_dtype = np.unicode

        actual_dtype = get_dtype_of_array(array)

        self.assertEqual(expected_dtype, actual_dtype)

    def test__get_dtype_of_array__unordered_integer_categories_return_as_expected(self):
        array = Series(data=[2, 3, 1, 3, 1, 2], dtype="category")
        expected_dtype = np.int32

        actual_dtype = get_dtype_of_array(array)

        self.assertEqual(expected_dtype, actual_dtype)

    def test__get_dtype_of_array__castable_dtypes_return_as_expected(self):
        types = [np.float64, np.int64]
        expected_dtypes = [np.float32, np.int32]

        for test_type_index in range(len(types)):
            with self.subTest(
                f"Testing get_dtype_of_array with castable type {types[test_type_index].__name__}", i=test_type_index
            ):
                array = Series(data=[], dtype=types[test_type_index])
                self.assertEqual(get_dtype_of_array(array), expected_dtypes[test_type_index])

    def test__get_dtype_of_array__unsupported_type_raises_exception(self):
        unsupported_array = Series(list([time() for _ in range(2)]), dtype="datetime64[ns]")

        with self.assertRaises(TypeError) as exception_context:
            get_dtype_of_array(unsupported_array)

        self.assertIn("unsupported", str(exception_context.exception))

    def test__get_schema_type_hint_of_array__supported_dtypes_return_as_expected(self):
        types = [np.float32, np.int32, np.bool_, str]
        expected_schema_hints = [{"type": "float32"}, {"type": "int32"}, {"type": "boolean"}, {"type": "string"}]

        for test_type_index in range(len(types)):
            with self.subTest(
                f"Testing get_schema_type_hint_of_array with type {types[test_type_index].__name__}", i=test_type_index
            ):
                array = Series(data=[], dtype=types[test_type_index])
                self.assertEqual(get_schema_type_hint_of_array(array), expected_schema_hints[test_type_index])

    def test__get_schema_type_hint_of_array__categories_return_as_expected(self):
        array = Series(data=["a", "b", "b"], dtype="category")
        expected_schema_hint = {"type": "categorical", "categories": ["a", "b"]}

        actual_schema_hint = get_schema_type_hint_of_array(array)

        self.assertEqual(expected_schema_hint, actual_schema_hint)

    def test__get_schema_type_hint_of_array__castable_dtypes_return_as_expected(self):
        types = [np.float64, np.int64]
        expected_schema_hints = [{"type": "float32"}, {"type": "int32"}]

        for test_type_index in range(len(types)):
            with self.subTest(
                f"Testing get_schema_type_hint_of_array with castable type {types[test_type_index].__name__}",
                i=test_type_index,
            ):
                array = Series(data=[], dtype=types[test_type_index])
                self.assertEqual(get_schema_type_hint_of_array(array), expected_schema_hints[test_type_index])

    def test__get_dtypes_and_schemas_of_dataframe__dtype_and_schema_returns_as_expected(self):
        float_array = Series(data=[1, 2, 3], dtype=np.dtype(np.float64))
        category_array = Series(data=["a", "b", "b"], dtype="category")
        dataframe = DataFrame({"float_array": float_array, "category_array": category_array})

        expected_data_types_dict = {"float_array": np.float32, "category_array": np.unicode}
        expected_schema_type_hints_dict = {
            "float_array": {"type": "float32"},
            "category_array": {"type": "categorical", "categories": ["a", "b"]},
        }

        actual_dataframe_data_types, actual_dataframe_schema_type_hints = get_dtypes_and_schemas_of_dataframe(dataframe)

        self.assertEqual(expected_data_types_dict, actual_dataframe_data_types)
        self.assertEqual(expected_schema_type_hints_dict, actual_dataframe_schema_type_hints)

    def test__convert_pandas_series_to_numpy__categorical_float64_to_float64_with_nans(self):
        expected_float_array = np.array([1.1, 2.2, np.NaN], dtype=np.float64)
        float_series = Series(data=[1.1, 2.2, np.NaN], dtype="category")

        actual_float_array = convert_pandas_series_to_numpy(float_series, np.float64)

        np.testing.assert_equal(expected_float_array, actual_float_array)

    def test__convert_pandas_series_to_numpy__float64_to_float64(self):
        expected_float_array = np.array([1.1, 2.2], dtype=np.float64)
        float_series = Series(data=[1.1, 2.2], dtype=np.dtype(np.float64))

        actual_float_array = convert_pandas_series_to_numpy(float_series, np.float64)

        np.testing.assert_equal(expected_float_array, actual_float_array)

    def test__convert_pandas_series_to_numpy__int64_to_int32_with_nans_throws_error(self):
        int_series = Series(data=[1, 2, np.NaN], dtype="category")

        with self.assertLogs(level="ERROR") as logger:
            convert_pandas_series_to_numpy(int_series, np.int32)

            self.assertIn(
                "Cannot convert a pandas Series object to an integer dtype if it contains NaNs", logger.output[0]
            )
