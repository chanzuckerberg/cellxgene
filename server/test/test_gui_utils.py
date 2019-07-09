import unittest

from server.gui.options_parser import parse_opt_string
from server.utils.errors import OptionsError

default_opts = {'title': None, 'layout': (), 'obs_names': None, 'var_names': None, 'max_category_items': 1000,
                'diffexp_lfc_cutoff': 0.01, 'label_file': None}


class UtilsTest(unittest.TestCase):

    def test_empty_opts(self):
        result = parse_opt_string("")
        self.assertEqual(result, default_opts)

    def test_one_opts(self):
        result = parse_opt_string("--title abcde")
        expected_opts = dict(default_opts)
        expected_opts["title"] = "abcde"
        self.assertEqual(result, expected_opts)

    def test_malformed_opts(self):
        with self.assertRaises(OptionsError):
            parse_opt_string("--ewwwww lkafsjkl")

    def test_complex_opts(self):
        result = parse_opt_string(
            "--title abcde --obs-names zzzzz --var-names 'xxx yyy' --max-category-items 234907 "
            "--diffexp-lfc-cutoff 0.00003847382")
        expected_opts = {'title': 'abcde', 'layout': (), 'obs_names': 'zzzzz', 'var_names': 'xxx yyy',
                         'max_category_items': 234907, 'diffexp_lfc_cutoff': 3.847382e-05,
                         'label_file': None}
        self.assertEqual(result, expected_opts)
