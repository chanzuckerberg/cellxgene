import unittest
import numpy as np

from server.common.utils.utils import (
    jsonify_strict,
)


class TestJsonifyStrict(unittest.TestCase):
    def test_jsonify_numpy_general_cases(self):
        self.assertEqual(jsonify_strict({}), "{}")
        self.assertEqual(jsonify_strict({"a": [], "b": "hello", "c": True}), '{"a": [], "b": "hello", "c": true}')

    def test_jsonify_numpy_float_edges(self):
        with self.assertRaises(ValueError):
            jsonify_strict({"nan": [np.nan]})

        with self.assertRaises(ValueError):
            jsonify_strict({"pinf": [np.inf]})

        with self.assertRaises(ValueError):
            jsonify_strict({"ninf": [np.inf]})

    def test_jsonify_numpy_ndarray(self):
        values = {
            "integer": [
                np.int8(0),
                np.int16(1),
                np.int32(2),
                np.int64(3),
                np.uint8(4),
                np.uint16(5),
                np.uint32(6),
                np.uint64(7),
            ],
            "floating": [
                np.float16(100.0),
                np.float32(101.0),
                np.float64(102.0),
            ],
        }
        # these just confirm our test assumptions
        self.assertTrue(isinstance(values["floating"][0], np.float16))
        self.assertTrue(isinstance(values["floating"][1], np.float32))
        self.assertTrue(isinstance(values["floating"][2], np.float64))
        self.assertTrue(isinstance(values["integer"][0], np.int8))
        self.assertTrue(isinstance(values["integer"][1], np.int16))
        self.assertTrue(isinstance(values["integer"][2], np.int32))
        self.assertTrue(isinstance(values["integer"][3], np.int64))
        self.assertTrue(isinstance(values["integer"][4], np.uint8))
        self.assertTrue(isinstance(values["integer"][5], np.uint16))
        self.assertTrue(isinstance(values["integer"][6], np.uint32))
        self.assertTrue(isinstance(values["integer"][7], np.uint64))
        # the actual test!
        self.assertEqual(
            jsonify_strict(values),
            '{"integer": [0, 1, 2, 3, 4, 5, 6, 7], "floating": [100.0, 101.0, 102.0]}',
        )
