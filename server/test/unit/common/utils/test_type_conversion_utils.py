import unittest
from unittest.mock import patch

import numpy as np
from pandas import Series

from server.common.utils.type_conversion_utils import can_cast_to_float32, can_cast_to_int32


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
