import json
from os import path
import unittest

from numpy import float32, int32
from werkzeug.datastructures import ImmutableMultiDict

from server.app.util.filter import _convert_variable, parse_filter, QueryStringError


class UtilTest(unittest.TestCase):
    """Test Case for endpoints"""

    def setUp(self):
        with open(path.join(path.dirname(__file__), "schema.json")) as fh:
            schema = json.load(fh)
        self.schema = schema["annotations"]

    def test_convert(self):
        five = _convert_variable("int32", "5")
        self.assertEqual(five, int32(5))

    def test_convert_zero(self):
        zero = _convert_variable("int32", "0")
        self.assertEqual(zero, 0)

    def test_convert_float(self):
        str_to_convert = "4.38719237129"
        val = _convert_variable("float32", str_to_convert)
        self.assertAlmostEqual(val, float32(str_to_convert))

    def test_convert_bool(self):
        str_to_convert = "false"
        val = _convert_variable("boolean", str_to_convert)
        self.assertFalse(val)
        str_to_convert = "true"
        val = _convert_variable("boolean", str_to_convert)
        self.assertTrue(val)
        str_to_convert = "0"
        with self.assertRaises(AssertionError):
            val = _convert_variable("boolean", str_to_convert)

    def test_empty_convert(self):
        empty = _convert_variable("int32", None)
        self.assertIsNone(empty)

    def test_bad_convert(self):
        with self.assertRaises(ValueError):
            _convert_variable("int32", "5.5")

    def test_bad_datatype(self):
        with self.assertRaises(AssertionError):
            _convert_variable("jkasdslkja", 1)

    def test_complex_filter(self):
        filter_dict = ImmutableMultiDict(
            [("obs:louvain", "NK cells"), ("obs:louvain", "CD8 T cells"), ("obs:n_counts", "3000,*")])
        filter_ = parse_filter(filter_dict, self.schema)
        self.assertIn("obs", filter_)
        self.assertEqual(filter_["obs"]["annotation_value"], [{"name": "louvain",
                                                               "values": ["NK cells", "CD8 T cells"]},
                                                              {"name": "n_counts",
                                                               "max": None, "min": 3000.0}])

    def test_bad_filter(self):
        bad_annotation_type = ImmutableMultiDict([("obs:tissue", "lung")])
        with self.assertRaises(QueryStringError):
            parse_filter(bad_annotation_type, self.schema)
        bad_axis = ImmutableMultiDict([("xyz:n_genes", "100,1000")])
        with self.assertRaises(QueryStringError):
            parse_filter(bad_axis, self.schema)

    def test_boolean_filter(self):
        schema = {
            "obs": [{"name": "bool_filter", "type": "boolean"}]
        }
        filter_dict = ImmutableMultiDict([("obs:bool_filter", "false")])
        filter_ = parse_filter(filter_dict, schema)
        self.assertIn("obs", filter_)
        self.assertEqual(filter_["obs"]["annotation_value"], [{"name": "bool_filter", "values": [False]}])
