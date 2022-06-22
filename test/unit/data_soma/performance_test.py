import timeit
import numpy as np

from test.unit import app_config
from test import PROJECT_ROOT, FIXTURES_ROOT
from test.unit.data_soma.dataset_handler import decompress_dataset
import test.decode_fbs as decode_fbs
from server.common.utils.data_locator import DataLocator
from server.data_soma.soma_adaptor import SomaAdaptor
from server.data_anndata.anndata_adaptor import AnndataAdaptor

if __name__ == "__main__":
    # run performance testing
    # best to run from root cellxgene dir and run `export PYTHONPATH=.` first

    NUMBER = 5
    soma_params = [
        (f"{PROJECT_ROOT}/example-dataset/tiledb-data/pbmc3k_processed.zip", "pbmc3k_processed", False, "auto"),
        # (f"{FIXTURES_ROOT}/tiledb-data/pbmc3k-CSC-gz_processed.zip", "pbmc3k-CSC-gz_processed", False, "auto"),
        # (f"{FIXTURES_ROOT}/tiledb-data/pbmc3k-CSR-gz_processed.zip", "pbmc3k-CSR-gz_processed", False, "auto"),
        # (f"{PROJECT_ROOT}/example-dataset/tiledb-data/pbmc3k_processed.zip", "pbmc3k_processed", True, "auto"),
        # (f"{FIXTURES_ROOT}/tiledb-data/pbmc3k-CSC-gz_processed.zip", "pbmc3k-CSC-gz_processed", True, "auto"),
        # (f"{FIXTURES_ROOT}/tiledb-data/pbmc3k-CSR-gz_processed.zip", "pbmc3k-CSR-gz_processed", True, "auto"),
        # (f"{PROJECT_ROOT}/example-dataset/tiledb-data/pbmc3k_processed.zip", "pbmc3k_processed", False, "normal"),
        # (f"{FIXTURES_ROOT}/tiledb-data/pbmc3k-CSC-gz_processed.zip", "pbmc3k-CSC-gz_processed", False, "normal"),
        # (f"{FIXTURES_ROOT}/tiledb-data/pbmc3k-CSR-gz_processed.zip", "pbmc3k-CSR-gz_processed", False, "normal"),
        # (f"{PROJECT_ROOT}/example-dataset/tiledb-data/pbmc3k_processed.zip", "pbmc3k_processed", True, "normal"),
        # (f"{FIXTURES_ROOT}/tiledb-data/pbmc3k-CSC-gz_processed.zip", "pbmc3k-CSC-gz_processed", True, "normal"),
        # (f"{FIXTURES_ROOT}/tiledb-data/pbmc3k-CSR-gz_processed.zip", "pbmc3k-CSR-gz_processed", True, "normal"),
        # (f"{FIXTURES_ROOT}/tiledb-data/pbmc3k_64_processed.zip", "pbmc3k_64_processed", False, "auto")
    ]
    anndata_params = [
        (f"{PROJECT_ROOT}/example-dataset/pbmc3k.h5ad", False, "auto"),
        # (f"{FIXTURES_ROOT}/pbmc3k-CSC-gz.h5ad", False, "auto"),
        # (f"{FIXTURES_ROOT}/pbmc3k-CSR-gz.h5ad", False, "auto"),
        # (f"{PROJECT_ROOT}/example-dataset/pbmc3k.h5ad", True, "auto"),
        # (f"{FIXTURES_ROOT}/pbmc3k-CSC-gz.h5ad", True, "auto"),
        # (f"{FIXTURES_ROOT}/pbmc3k-CSR-gz.h5ad", True, "auto"),
        # (f"{PROJECT_ROOT}/example-dataset/pbmc3k.h5ad", False, "normal"),
        # (f"{FIXTURES_ROOT}/pbmc3k-CSC-gz.h5ad", False, "normal"),
        # (f"{FIXTURES_ROOT}/pbmc3k-CSR-gz.h5ad", False, "normal"),
        # (f"{PROJECT_ROOT}/example-dataset/pbmc3k.h5ad", True, "normal"),
        # (f"{FIXTURES_ROOT}/pbmc3k-CSC-gz.h5ad", True, "normal"),
        # (f"{FIXTURES_ROOT}/pbmc3k-CSR-gz.h5ad", True, "normal"),
        # (f"{FIXTURES_ROOT}/pbmc3k_64.h5ad", False, "auto"),
    ]

    soma_times = []
    anndata_times = []
    for soma_param, anndata_param in zip(soma_params, anndata_params):
        soma_path, soma_name, soma_backed, soma_x = soma_param
        anndata_path, anndata_backed, anndata_x = anndata_param

        soma_path = decompress_dataset(soma_path, soma_name)
        soma_config = app_config(
            soma_path,
            soma_backed,
            extra_dataset_config=dict(X_approximate_distribution=soma_x),
        )
        anndata_config = app_config(
            anndata_path,
            anndata_backed,
            extra_dataset_config=dict(X_approximate_distribution=anndata_x),
        )
        data_soma = SomaAdaptor(DataLocator(soma_path), soma_config)
        data_anndata = AnndataAdaptor(DataLocator(anndata_path), anndata_config)

        time_soma_load = timeit.timeit("SomaAdaptor(DataLocator(soma_path), soma_config)", 'from __main__ import SomaAdaptor, DataLocator, soma_path, soma_config', number=NUMBER) / NUMBER
        time_anndata_load = timeit.timeit("AnndataAdaptor(DataLocator(anndata_path), anndata_config)", 'from __main__ import AnndataAdaptor, DataLocator, anndata_path, anndata_config', number=NUMBER) / NUMBER

        time_soma_embednames = timeit.timeit("data_soma.get_embedding_names()", 'from __main__ import data_soma', number=NUMBER) / NUMBER
        time_anndata_embednames = timeit.timeit("data_anndata.get_embedding_names()", 'from __main__ import data_anndata', number=NUMBER) / NUMBER

        time_soma_embedarr = timeit.timeit("data_soma.get_embedding_array('pca')", 'from __main__ import data_soma', number=NUMBER) / NUMBER
        time_anndata_embedarr = timeit.timeit("data_anndata.get_embedding_array('pca')", 'from __main__ import data_anndata', number=NUMBER) / NUMBER

        filter_ = {
            "filter": {"var": {"annotation_value": [{"name": "n_cells", "min": 10}], "index": [1, 99, [200, 300]]}}
        }
        time_soma_filter = timeit.timeit("data_soma.data_frame_to_fbs_matrix(filter_['filter'], 'var')", 'from __main__ import data_soma, filter_', number=NUMBER) / NUMBER
        time_anndata_filter = timeit.timeit("data_anndata.data_frame_to_fbs_matrix(filter_['filter'], 'var')", 'from __main__ import data_anndata, filter_', number=NUMBER) / NUMBER

        time_soma_xdist = timeit.timeit("data_soma.get_X_approximate_distribution()", 'from __main__ import data_soma', number=NUMBER) / NUMBER
        time_anndata_xdist = timeit.timeit("data_anndata.get_X_approximate_distribution()", 'from __main__ import data_anndata', number=NUMBER) / NUMBER


        coll_soma = [time_soma_load, time_soma_embednames, time_soma_embedarr, time_soma_filter, time_soma_xdist]
        coll_anndata = [time_anndata_load, time_anndata_embednames, time_anndata_embedarr, time_anndata_filter, time_anndata_xdist]
        soma_times.append(coll_soma)
        anndata_times.append(coll_anndata)


    print("\nSOMA vs Anndata PERFORMANCE TESTS")
    print("=================")

    for soma, anndata, soma_params, anndata_params in zip(soma_times, anndata_times, soma_params, anndata_params):
        print(f"SOMA    - {soma}")
        print(f"Anndata - {anndata}")
        print(f"Tests - [load, embednames, embedarr, filter, xdist]")
        print(f"SOMA params = {soma_params}, Anndata params = {anndata_params}")
        print("----------")
