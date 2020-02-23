import shutil
import tempfile
from os import path

import pandas as pd

from server.common.annotations import AnnotationsLocalFile
from server.common.data_locator import DataLocator
from server.data_common.fbs.matrix import encode_matrix_fbs
from server.data_common.matrix_loader import MatrixDataLoader, MatrixDataType


def data_with_tmp_annotations(ext: MatrixDataType, annotations_fixture=False):
    tmp_dir = tempfile.mkdtemp()
    annotations_file = path.join(tmp_dir, "test_annotations.csv")
    if annotations_fixture:
        shutil.copyfile(f"test/test_datasets/pbmc3k-annotations.csv", annotations_file)
    args = {
        "layout": ["umap"],
        "max_category_items": 100,
        "obs_names": None,
        "var_names": None,
        "diffexp_lfc_cutoff": 0.01,
    }
    fname = {
        MatrixDataType.H5AD: "../example-dataset/pbmc3k.h5ad",
        MatrixDataType.CXG: "test/test_datasets/pbmc3k.cxg",
    }[ext]
    data_locator = DataLocator(fname)
    data = MatrixDataLoader(data_locator.abspath()).open(args)
    annotations = AnnotationsLocalFile(None, annotations_file)
    return data, tmp_dir, annotations


def make_fbs(data):
    df = pd.DataFrame(data)
    return encode_matrix_fbs(matrix=df, row_idx=None, col_idx=df.columns)


def skip_if(condition, reason: str):
    def decorator(f):
        def wraps(self, *args, **kwargs):
            if condition(self):
                self.skipTest(reason)
            else:
                f(self, *args, **kwargs)

        return wraps

    return decorator
