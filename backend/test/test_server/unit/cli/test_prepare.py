import unittest

import pandas as pd

from backend.server.cli.prepare import make_index_unique


class CLIPrepareTests(unittest.TestCase):
    """ Test cases for CLI prepare logic """

    def test_make_index_unique(self):
        index = pd.Index(["SNORD113", "SNORD113", "SNORD113-1"])
        result = make_index_unique(index)
        expected = pd.Index(["SNORD113", "SNORD113-2", "SNORD113-1"])
        self.assertTrue(all(left == right for left, right in zip(result.values, expected.values)))
