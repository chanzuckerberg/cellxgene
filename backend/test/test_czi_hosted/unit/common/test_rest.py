import unittest
from urllib.parse import parse_qs
from werkzeug.datastructures import MultiDict
from backend.common.errors import FilterError
from backend.czi_hosted.common.rest import _query_parameter_to_filter


def _qsparse(qs):
    """emulate what Flask/Werkzeug do to our QS"""
    return MultiDict(parse_qs(qs))


class FilterParseTests(unittest.TestCase):
    """Test cases for various filter parsing"""

    def test_queryparam_to_filter_parse(self):
        # categories
        self.assertEqual(
            _query_parameter_to_filter(_qsparse("obs:foo=bar&var:baz=133&var:baz=A&obs:baz=foo")),
            {
                "obs": {"annotation_value": [{"name": "foo", "values": ["bar"]}, {"name": "baz", "values": ["foo"]}]},
                "var": {"annotation_value": [{"name": "baz", "values": ["133", "A"]}]},
            },
        )

        # ranges
        self.assertEqual(
            _query_parameter_to_filter(_qsparse("obs:A=1,99&obs:B=*,100&obs:C=0,*")),
            {
                "obs": {
                    "annotation_value": [
                        {"name": "A", "min": 1, "max": 99.0},
                        {"name": "B", "max": 100.0},
                        {"name": "C", "min": 0.0},
                    ]
                },
            },
        )

        # combo
        self.assertEqual(
            _query_parameter_to_filter(_qsparse("var:B=YES&var:A=1,99&var:B=NO")),
            {
                "var": {
                    "annotation_value": [
                        {"name": "B", "values": ["YES", "NO"]},
                        {"name": "A", "min": 1.0, "max": 99.0},
                    ]
                },
            },
        )

    def test_queryparam_to_filter_escaping(self):
        self.assertEqual(
            _query_parameter_to_filter(_qsparse("obs:var=%2521%252C%253AOK%253D&obs:A%2521=YO")),
            {"obs": {"annotation_value": [{"name": "var", "values": ["!,:OK="]}, {"name": "A!", "values": ["YO"]}]}},
        )

    def test_queryparam_to_filter_errors(self):

        # should raise FilterError
        filter_errors = [
            "foo=bar",  # no axis
            "X=&Y=3",  # no value
            "X&Y=3",  # no value
            "moo:foo=bar",  # bad axis
            "obs:x=1,A",  # non-numeric range
            "var:X=1,2&var:X=3,4",  # duplicate ranges
            "var:Y=,",
            "var:Y=2,",
            "var:Y=,5",
            "var:Y=*,",
            "var:Y=,*",
            "var:Y=*,*",
        ]

        for qs in filter_errors:
            with self.assertRaises(FilterError):
                _query_parameter_to_filter(_qsparse(qs))
