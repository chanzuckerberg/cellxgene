import zipfile
import os
from test import PROJECT_ROOT

def decompress_dataset(path: str, name: str):
    """Decompresses SOMA datasets for testing use"""
    # first check if dataset with name exists
    if os.path.isdir(f"{PROJECT_ROOT}/test/fixtures/tiledb-data/decompressed/{name}"):
        return f"{PROJECT_ROOT}/test/fixtures/tiledb-data/decompressed/{name}"

    if zipfile.is_zipfile(path):
        # decompress into new directory and return new path
        with zipfile.ZipFile(path, 'r') as zip_ref:
            zip_ref.extractall(f"{PROJECT_ROOT}/test/fixtures/tiledb-data/decompressed")
        return f"{PROJECT_ROOT}/test/fixtures/tiledb-data/decompressed/{name}"
    else:
        return path