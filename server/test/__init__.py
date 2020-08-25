import os
import random
import shutil
import string
import tempfile
import time
from contextlib import contextmanager
from os import path, popen
from subprocess import Popen

import pandas as pd
import requests

from server.common.annotations.hosted_tiledb import AnnotationsHostedTileDB
from server.common.annotations.local_file_csv import AnnotationsLocalFile
from server.common.app_config import AppConfig, DEFAULT_SERVER_PORT
from server.common.data_locator import DataLocator
from server.common.utils.utils import find_available_port
from server.data_common.fbs.matrix import encode_matrix_fbs
from server.data_common.matrix_loader import MatrixDataLoader, MatrixDataType
from server.db.db_utils import DbUtils

PROJECT_ROOT = popen("git rev-parse --show-toplevel").read().strip()
FIXTURES_ROOT = PROJECT_ROOT + "/server/test/fixtures"


def data_with_tmp_tiledb_annotations(ext: MatrixDataType):
    tmp_dir = tempfile.mkdtemp()
    fname = {
        MatrixDataType.H5AD: f"{PROJECT_ROOT}/example-dataset/pbmc3k.h5ad",
        MatrixDataType.CXG: "test/fixtures/pbmc3k.cxg",
    }[ext]
    data_locator = DataLocator(fname)
    config = AppConfig()
    config.update_server_config(
        multi_dataset__dataroot=data_locator.path,
        authentication__type="test",
    )
    config.update_default_dataset_config(
        embeddings__names=["umap"],
        presentation__max_categories=100,
        diffexp__lfc_cutoff=0.01,
        user_annotations__type="hosted_tiledb_array",
        user_annotations__hosted_tiledb_array__db_uri="postgresql://postgres:test_pw@localhost:5432",
        user_annotations__hosted_tiledb_array__hosted_file_directory=tmp_dir
    )

    config.complete_config()

    data = MatrixDataLoader(data_locator.abspath()).open(config)
    annotations = AnnotationsHostedTileDB(
        tmp_dir,
        DbUtils("postgresql://postgres:test_pw@localhost:5432"),
    )
    return data, tmp_dir, annotations


def data_with_tmp_annotations(ext: MatrixDataType, annotations_fixture=False):
    tmp_dir = tempfile.mkdtemp()
    annotations_file = path.join(tmp_dir, "test_annotations.csv")
    if annotations_fixture:
        shutil.copyfile(f"{PROJECT_ROOT}/server/test/fixtures/pbmc3k-annotations.csv", annotations_file)
    fname = {
        MatrixDataType.H5AD: f"{PROJECT_ROOT}/example-dataset/pbmc3k.h5ad",
        MatrixDataType.CXG: "test/fixtures/pbmc3k.cxg",
    }[ext]
    data_locator = DataLocator(fname)
    config = AppConfig()
    config.update_server_config(
        single_dataset__obs_names=None, single_dataset__var_names=None, single_dataset__datapath=data_locator.path
    )
    config.update_default_dataset_config(
        embeddings__names=["umap"], presentation__max_categories=100, diffexp__lfc_cutoff=0.01,
    )

    config.complete_config()
    data = MatrixDataLoader(data_locator.abspath()).open(config)
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


def app_config(data_locator, backed=False, extra_server_config={}, extra_dataset_config={}):
    config = AppConfig()
    config.update_server_config(
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


def random_string(n):
    return "".join(random.choice(string.ascii_letters) for _ in range(n))


def start_test_server(command_line_args=[], app_config=None):
    """
    Command line arguments can be passed in, as well as an app_config.
    This function is meant to be used like this, for example:

    with test_server(...) as server:
        r = requests.get(f"{server}/...")
        // check r

    where the server can be accessed within the context, and is terminated when
    the context is exited.
    The port is automatically set using find_available_port.
    The verbose flag is automatically set to True.
    If an app_config is provided, then this function writes a temporary
    yaml config file, which this server will read and parse.
    """

    start = random.randint(DEFAULT_SERVER_PORT, 2 ** 16 - 1)
    port = int(os.environ.get("CXG_SERVER_PORT", start))
    port = find_available_port("localhost", port)
    command = ["cellxgene", "--no-upgrade-check", "launch", "--verbose", "--port=%d" % port] + command_line_args

    tempdir = None
    if app_config:
        tempdir = tempfile.TemporaryDirectory()
        config_file = os.path.join(tempdir.name, "config.yaml")
        app_config.write_config(config_file)
        command.extend(["-c", config_file])

    server = f"http://localhost:{port}"
    ps = Popen(command)

    for _ in range(10):
        try:
            requests.get(f"{server}/health")
            break
        except requests.exceptions.ConnectionError:
            time.sleep(1)

    if tempdir:
        tempdir.cleanup()

    return ps, server


def stop_test_server(ps):
    try:
        ps.terminate()
    except ProcessLookupError:
        pass


@contextmanager
def test_server(command_line_args=[], app_config=None):
    """A context to run the cellxgene server."""

    ps, server = start_test_server(command_line_args, app_config)
    try:
        yield server
    finally:
        try:
            stop_test_server(ps)
        except ProcessLookupError:
            pass
