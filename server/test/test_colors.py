import unittest

import anndata
from server.common.colors import convert_color_to_hex_format, convert_anndata_category_colors_to_cxg_category_colors
from server.test import PROJECT_ROOT
from server.test.test_datasets.fixtures import pbmc3k_colors


class ColorsTest(unittest.TestCase):
    """ Test color helper functions """

    def test_convert_color_to_hex_format(self):
        self.assertEqual(convert_color_to_hex_format("wheat"), "#f5deb3")
        self.assertEqual(convert_color_to_hex_format((245, 222, 179)), "#f5deb3")
        self.assertEqual(convert_color_to_hex_format([245, 222, 179]), "#f5deb3")
        self.assertEqual(convert_color_to_hex_format("#f5deb3"), "#f5deb3")
        self.assertEqual(
            convert_color_to_hex_format([0.9607843137254902, 0.8705882352941177, 0.7019607843137254]), "#f5deb3"
        )

    def test_anndata_colors_to_cxg_colors(self):
        adata = anndata.read_h5ad(f"{PROJECT_ROOT}/example-dataset/pbmc3k.h5ad")
        self.assertEqual(convert_anndata_category_colors_to_cxg_category_colors(adata), pbmc3k_colors)
