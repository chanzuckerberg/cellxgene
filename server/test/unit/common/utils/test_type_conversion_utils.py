import unittest
from unittest.mock import patch

import numpy as np
from pandas import Series

from server.common.utils.type_conversion_utils import can_cast_to_float32, can_cast_to_int32, get_dtype_of_array


class TestTypeConversionUtils(unittest.TestCase):

    def test__can_cast_to_float32__string_is_false(self):
        array_to_convert = Series(data=["1", "2", "3"], dtype=str)

        can_cast = can_cast_to_float32(array_to_convert)

        self.assertFalse(can_cast)

    def test__can_cast_to_float32__int_is_true_warning_outputted(self):
        array_to_convert = Series(data=[1, 2, 3], dtype=np.dtype(np.float64))

        with self.assertLogs(level="WARN") as logger:
            can_cast = can_cast_to_float32(array_to_convert)
            self.assertIn(f"may lose precision", logger.output[0])

        self.assertTrue(can_cast)

    @patch("logging.warning")
    def test__can_cast_to_float64__int_is_false(self, mock_log_warning):
        array_to_convert = Series(data=[1, 2, 3], dtype=np.dtype(np.float32))

        can_cast = can_cast_to_float32(array_to_convert)

        self.assertTrue(can_cast)
        assert not mock_log_warning.called

    def test__can_cast_to_int32__string_is_false(self):
        array_to_convert = Series(data=["1", "2", "3"], dtype=str)

        can_cast = can_cast_to_int32(array_to_convert)

        self.assertFalse(can_cast)

    def test__can_cast_to_int32__int64_is_true(self):
        array_to_convert = Series(data=["1", "2", "3"], dtype=np.dtype(np.int64))

        can_cast = can_cast_to_int32(array_to_convert)

        self.assertTrue(can_cast)

    def test__can_cast_to_int32__int16_is_true(self):
        array_to_convert = Series(data=["1", "2", "3"], dtype=np.dtype(np.int16))

        can_cast = can_cast_to_int32(array_to_convert)

        self.assertTrue(can_cast)

    def test__can_cast_to_int32__int64_with_large_value_is_false(self):
        array_to_convert = Series(data=["3000000000", "2", "3"], dtype=np.dtype(np.int64))

        can_cast = can_cast_to_int32(array_to_convert)

        self.assertFalse(can_cast)

    def test__get_dtype_of_array__supported_dtypes_return_as_expected(self):
        types = [np.float32, np.int32, np.bool_, np.str]
        expected_dtypes = [np.float32, np.int32, np.uint8, np.unicode]

        for test_type_index in range(len(types)):
            with self.subTest(f"Testing get_dtype_of_array with type {types[test_type_index].__name__}",
                              i=test_type_index):
                array = Series(data=[], dtype=np.dtype(types[test_type_index]))
                self.assertEqual(get_dtype_of_array(array), expected_dtypes[test_type_index])
