import unittest

from unittest.mock import MagicMock
import sys
sys.path.insert(0, "../app")
from util.filter import _convert_variable, parse_filter


class UtilTest(unittest.TestCase):
    """Test Case for endpoints"""

    def setUp(self):
       self.schema = {
            "cluster": {
                "displayname": "Cluster",
                "include": True,
                "type": "int",
                "variabletype": "categorical"
            },
            "louvain": {
                "displayname": "Louvain Cluster",
                "include": True,
                "type": "string",
                "variabletype": "categorical"
            },
            "n_genes": {
                "displayname": "Num Genes",
                "include": True,
                "type": "int",
                "variabletype": "continuous"
            }
        }

    def test_convert(self):
        five = _convert_variable("int", "5")
        assert five == 5

    def test_convert_zero(self):
        zero = _convert_variable("int", "0")
        assert zero == 0

    def test_empty_convert(self):
        empty = _convert_variable("int", None)
        assert empty is None

    def test_bad_convert(self):
        with self.assertRaises(ValueError):
            _convert_variable("int", "5.5")

    def test_filter_categorical(self):
        filterMock = MagicMock()
        filterMock.__iter__.return_value = iter(["louvain"])
        filterMock.getlist.return_value = ["B cells", "T cells"]
        query = parse_filter(filterMock, self.schema)
        assert query == {"louvain": {"variable_type": "categorical", "value_type": "string", "query": ["B cells", "T cells"]}}
        filterMock.__iter__.return_value = iter(["cluster"])
        filterMock.getlist.return_value = ["1", "2"]
        query = parse_filter(filterMock, self.schema)
        assert query == {"cluster": {"variable_type": "categorical", "value_type": "int", "query": [1, 2]}}

    def test_filter_contiunous(self):
        filterMock = MagicMock()
        filterMock.__iter__.return_value = iter(["n_genes"])
        filterMock.getlist.return_value = ["0,100"]
        query = parse_filter(filterMock, self.schema)
        assert query == {"n_genes": {"variable_type": "continuous", "value_type": "int", "query": {"min": 0, "max": 100}}}
        filterMock.__iter__.return_value = iter(["n_genes"])
        filterMock.getlist.return_value = ["*,100"]
        query = parse_filter(filterMock, self.schema)
        assert query == {"n_genes": {"variable_type": "continuous", "value_type": "int", "query": {"min": None, "max": 100}}}
        filterMock.__iter__.return_value = iter(["n_genes"])
        filterMock.getlist.return_value = ["0,*"]
        query = parse_filter(filterMock, self.schema)
        assert query == {"n_genes": {"variable_type": "continuous", "value_type": "int", "query": {"min": 0, "max": None}}}
