import json
import unittest

import numpy as np
from werkzeug.datastructures import MultiDict

from server.common.rest import _query_parameter_to_filter
from server.common.utils.data_locator import DataLocator
from server.data_cxg.cxg_dataset import CxgDataset
from test import FIXTURES_ROOT, decode_fbs
from test.fixtures.fixtures import pbmc3k_colors
from test.unit import app_config


class TestCxgDataset(unittest.TestCase):
    def test_get_colors(self):
        data = self.get_data("pbmc3k.cxg")
        self.assertDictEqual(data.get_colors(), pbmc3k_colors)
        data = self.get_data("pbmc3k_v0.cxg")
        self.assertDictEqual(data.get_colors(), dict())

    def get_data(self, fixture):
        data_locator = f"{FIXTURES_ROOT}/{fixture}"
        config = app_config(data_locator)
        return CxgDataset(DataLocator(data_locator), config)

    def test_schema_filters_int64(self):
        data = self.get_data("dataset_with_int64.cxg")
        schema = data.get_schema()
        with open(f"{FIXTURES_ROOT}/schema_without_int64.json", "r") as json_file:
            expected_schema = json.load(json_file)
            self.assertEqual(json.dumps(schema), json.dumps(expected_schema))

    def test_tdb_bug(self):
        """
        This gives different results on 0.12.4 vs 0.13.1. Reported to TileDB
        and fixed in 0.13.2. Test case remains in case of regression.
        """
        for dataset in ["pbmc3k.cxg", "pbmc3k_sparse.cxg"]:
            with self.subTest(dataset=dataset):
                data = self.get_data(dataset)
                filt = _query_parameter_to_filter(
                    MultiDict(
                        [
                            ("var:name_0", "F5"),
                            ("var:name_0", "BEB3"),
                            ("var:name_0", "SIK1"),
                        ]
                    )
                )
                dat = data.summarize_var("mean", filt, 0)
                summary = decode_fbs.decode_matrix_FBS(dat)
                self.assertDictContainsSubset({"n_rows": 2638, "n_cols": 1, "row_idx": None}, summary)
                self.assertIs(type(summary["columns"]), list)
                self.assertEqual(len(summary["columns"]), 1)
                self.assertEqual(len(summary["columns"][0]), 2638)
                self.assertEqual(summary["columns"][0].sum(), np.float32(-19.00301))

    def test_tdb_bug_lossy(self):
        """
        This gives different results on 0.12.4 vs 0.13.1. Reported to TileDB
        and fixed in 0.13.2. Test case remains in case of regression.
        """
        data = self.get_data("pbmc3k.cxg")
        filt = _query_parameter_to_filter(
            MultiDict(
                [
                    ("var:name_0", "F5"),
                    ("var:name_0", "BEB3"),
                    ("var:name_0", "SIK1"),
                ]
            )
        )
        dat = data.summarize_var("mean", filt, 0, 500)
        summary = decode_fbs.decode_matrix_FBS(dat)
        self.assertDictContainsSubset({"n_rows": 2638, "n_cols": 1, "row_idx": None}, summary)
        self.assertIs(type(summary["columns"]), list)
        self.assertEqual(len(summary["columns"]), 1)
        self.assertEqual(len(summary["columns"][0]), 2638)
        self.assertEqual(summary["columns"][0].sum(), np.float32(-35.90116))
