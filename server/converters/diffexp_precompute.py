"""This program precomputes data to speed up diffexp."""

import sys
import os
import tiledb
import argparse
import numpy as np
import shutil
from server.common.app_config import AppConfig
from server.data_common.matrix_loader import MatrixDataLoader
from server.compute.diffexp_generic import compute_indices_hash, compute_indices_hash_version
from server.compute.diffexp_cxg import compute_mean_var
from server.data_cxg.cxg_adaptor import CxgAdaptor
from server.common.errors import DatasetAccessError, ComputeError


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("dataset", help="input CXG dataset")
    parser.add_argument("--max_groups", "--max", type=int, default=100)
    parser.add_argument("--min_groups", "--min", type=int, default=2)
    parser.add_argument("--max_rows", "--maxr", type=int, default=1000000)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    dataset = args.dataset

    try:
        app_config = AppConfig()
        app_config.single_dataset__datapath = dataset
        app_config.server__verbose = True
        app_config.diffexp__alg_cxg__target_workunit = 8000000
        app_config.complete_config()

        loader = MatrixDataLoader(dataset)
        adaptor = loader.open(app_config)
    except DatasetAccessError as e:
        print(f"Unable to open {dataset}: {str(e)}")
        sys.exit(1)

    if not isinstance(adaptor, CxgAdaptor):
        print(f"Expected {dataset} to be a valid CxgAdaptor")
        sys.exit(1)

    if not adaptor.data_locator.islocal():
        print(f"Expected {dataset} to be local")
        sys.exit(1)

    diffexp_group = os.path.join(dataset, "diffexp")
    if os.path.exists(diffexp_group):
        print(f"{diffexp_group} exists")
        if args.overwrite:
            print(f"{diffexp_group} exists, will overwrite")
            shutil.rmtree(diffexp_group)
        else:
            print(f"{diffexp_group} exists, use '--overwrite' to overwrite existing files")
            sys.exit(1)

    shape = adaptor.get_shape()
    attrs = adaptor.get_obs_columns()

    # determine which obs values to precompute, and how many.
    num_precompute = 0
    attr_precompute = []
    for attr in attrs:
        values = adaptor.query_obs_array(attr)
        uvalues = np.unique(values)
        ngroups = len(uvalues)
        if ngroups > args.max_groups:
            print(f"Skipping '{attr}': too many groups {ngroups} > {args.max_groups}")
            continue
        if ngroups < args.min_groups:
            print(f"Skipping '{attr}': too few groups {ngroups} < {args.min_groups}")
            continue

        print(f"Queueing '{attr}': num groups {ngroups}")
        num_precompute += ngroups
        attr_precompute.append(attr)

    hashkey = np.empty((num_precompute,), dtype=object)
    metakey = np.empty((num_precompute, 2), dtype=object)
    mean_var = np.zeros((2, num_precompute, shape[1]), dtype=np.float64)

    print()
    num_precompute = precompute(adaptor, attr_precompute, hashkey, metakey, mean_var, args.max_rows)

    np.resize(hashkey, (num_precompute,))
    np.resize(metakey, (num_precompute, 2))
    np.resize(mean_var, (2, num_precompute, shape[1]))

    # check if group exists
    tiledb.group_create(diffexp_group, CxgAdaptor.tiledb_ctx)
    hash_version = compute_indices_hash_version()
    save_array(diffexp_group, "hashkey", hashkey, dict(hash_version=hash_version))
    save_array(diffexp_group, "metakey", metakey)
    save_array(diffexp_group, "mean_var", mean_var)

    # TODO save metadata about the hash function version, in case it gets changed


def precompute(adaptor, attrnames, hashkey, metakey, mean_var, max_rows):
    hash_version = compute_indices_hash_version()
    pindex = 0
    for attr in attrnames:
        values = adaptor.query_obs_array(attr)
        uvalues = np.unique(values)
        ngroups = len(uvalues)
        print(f"Processing '{attr}': num groups = {ngroups}")
        for i in range(ngroups):
            val = uvalues[i]
            groupindices = np.where(values == val)[0]
            grouplen = len(groupindices)
            if grouplen > max_rows:
                print(f" '{attr}'='{str(val)}': num indices = {grouplen} > {max_rows}: skipping")
                continue
            else:
                print(f" '{attr}'='{str(val)}': num indices = {grouplen}")
            try:
                tmean, tvariance = compute_mean_var(adaptor, groupindices)
            except ComputeError as e:
                print(f"   error: {str(e)}")
                continue
            mean_var[0][pindex] = tmean
            mean_var[1][pindex] = tvariance
            hashkey[pindex] = compute_indices_hash(groupindices, hash_version)
            metakey[pindex][0] = attr
            metakey[pindex][1] = str(val)
            pindex += 1

    return pindex


def save_array(diffexp_group, name, array, meta=None):
    varname = os.path.join(diffexp_group, name)
    A = tiledb.DenseArray.from_numpy(varname, array, ctx=CxgAdaptor.tiledb_ctx)
    if meta:
        with tiledb.DenseArray(varname, mode="w") as A:
            for k, v in meta.items():
                A.meta[k] = v

    print(f"Write {varname}: shape={array.shape} dtype={array.dtype}")


if __name__ == "__main__":
    main()
