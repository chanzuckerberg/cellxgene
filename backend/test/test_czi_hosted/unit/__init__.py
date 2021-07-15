import logging
import shutil
import tempfile
import unittest

from os import path

import pandas as pd
from flask_compress import Compress
from flask_cors import CORS

from backend.czi_hosted.common.annotations.hosted_tiledb import AnnotationsHostedTileDB
from backend.czi_hosted.common.annotations.local_file_csv import AnnotationsLocalFile
from backend.czi_hosted.common.config.app_config import AppConfig
from backend.common.utils.data_locator import DataLocator
from backend.common.fbs.matrix import encode_matrix_fbs
from backend.czi_hosted.data_common.matrix_loader import MatrixDataType, MatrixDataLoader
from backend.czi_hosted.db.db_utils import DbUtils
from backend.czi_hosted.app.app import Server
from backend.test import PROJECT_ROOT, FIXTURES_ROOT


def data_with_tmp_tiledb_annotations(ext: MatrixDataType):
    tmp_dir = tempfile.mkdtemp()
    fname = {
        MatrixDataType.H5AD: f"{PROJECT_ROOT}/example-dataset/pbmc3k.h5ad",
        MatrixDataType.CXG: "test/fixtures/pbmc3k.cxg",
    }[ext]
    data_locator = DataLocator(fname)
    config = AppConfig()
    config.update_server_config(
        app__flask_secret_key="secret",
        multi_dataset__dataroot=data_locator.path,
        authentication__type="test",
        authentication__insecure_test_environment=True,
    )
    config.update_default_dataset_config(
        embeddings__names=["umap"],
        presentation__max_categories=100,
        diffexp__lfc_cutoff=0.01,
        user_annotations__type="hosted_tiledb_array",
        user_annotations__hosted_tiledb_array__db_uri="postgresql://postgres:test_pw@localhost:5432",
        user_annotations__hosted_tiledb_array__hosted_file_directory=tmp_dir,
    )

    config.complete_config()

    data = MatrixDataLoader(data_locator.abspath()).open(config)
    annotations = AnnotationsHostedTileDB(
        {
            "user-annotations": True,
            "genesets-save": False,
        },
        tmp_dir,
        DbUtils("postgresql://postgres:test_pw@localhost:5432"),
    )
    return data, tmp_dir, annotations


def data_with_tmp_annotations(ext: MatrixDataType, annotations_fixture=False):
    tmp_dir = tempfile.mkdtemp()
    annotations_file = path.join(tmp_dir, "test_annotations.csv")
    if annotations_fixture:
        shutil.copyfile(f"{FIXTURES_ROOT}/pbmc3k-annotations.csv", annotations_file)
    fname = {
        MatrixDataType.H5AD: f"{PROJECT_ROOT}/example-dataset/pbmc3k.h5ad",
        MatrixDataType.CXG: f"{FIXTURES_ROOT}/pbmc3k.cxg",
    }[ext]
    data_locator = DataLocator(fname)
    config = AppConfig()
    config.update_server_config(
        app__flask_secret_key="secret",
        single_dataset__obs_names=None,
        single_dataset__var_names=None,
        single_dataset__datapath=data_locator.path,
    )
    config.update_default_dataset_config(
        embeddings__names=["umap"],
        presentation__max_categories=100,
        diffexp__lfc_cutoff=0.01,
    )

    config.complete_config()
    data = MatrixDataLoader(data_locator.abspath()).open(config)
    annotations = AnnotationsLocalFile(
        {
            "user-annotations": True,
            "genesets-save": False,
        },
        None,
        annotations_file,
    )
    return data, tmp_dir, annotations, config


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


def app_config(data_locator, backed=False, extra_server_config={}, extra_dataset_config={}):
    config = AppConfig()
    config.update_server_config(
        app__flask_secret_key="secret",
        single_dataset__obs_names=None,
        single_dataset__var_names=None,
        adaptor__anndata_adaptor__backed=backed,
        single_dataset__datapath=data_locator,
        limits__diffexp_cellcount_max=None,
        limits__column_request_max=None,
    )
    config.update_default_dataset_config(
        embeddings__names=["umap", "tsne", "pca"], presentation__max_categories=100, diffexp__lfc_cutoff=0.01
    )
    config.update_server_config(**extra_server_config)
    config.update_default_dataset_config(**extra_dataset_config)
    config.complete_config()
    return config


class TestServer(Server):
    def __init__(self, app_config):
        super().__init__(app_config)

    @staticmethod
    def _before_adding_routes(app, app_config):
        app.config["COMPRESS_MIMETYPES"] = [
            "text/html",
            "text/css",
            "text/xml",
            "application/json",
            "application/javascript",
            "application/octet-stream",
        ]
        Compress(app)
        if app_config.server_config.app__debug:
            CORS(app, supports_credentials=True)


class BaseTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls, app_config=None):
        cls.TEST_URL_BASE = "/d/pbmc3k.cxg/api/v0.2/"
        cls.maxDiff = None
        cls.app = cls.create_app(app_config)

    @classmethod
    def create_app(cls, app_config=None):
        if not app_config:
            app_config = AppConfig()
            app_config.update_server_config(
                authentication__type="test",
                authentication__insecure_test_environment=True,
                app__flask_secret_key="testing",
                app__debug=True,
                multi_dataset__dataroot=f"{FIXTURES_ROOT}",
                multi_dataset__index=True,
                multi_dataset__allowed_matrix_types=["cxg"]
            )
        app_config.complete_config(logging.info)

        app = TestServer(app_config).app

        app.testing = True
        app.debug = True

        return app
