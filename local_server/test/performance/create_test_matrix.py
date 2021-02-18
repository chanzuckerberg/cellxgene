import anndata
import argparse
import random
import scipy
import numpy as np


def main():
    parser = argparse.ArgumentParser("A command to generate test h5ad files")
    parser.add_argument("output", help="Name of the output file")
    parser.add_argument("nobs", type=int, help="Number of observations (rows)")
    parser.add_argument("nvar", type=int, help="Number of variables (columns)")
    parser.add_argument("-n", "--nnz-percent", type=float, default=100, help="percent of non-zeros")
    parser.add_argument("-c", "--col-shift", action="store_true", help="add a random value to each column")
    parser.add_argument("--seed", type=int, default=None, help="add a random value to each column")

    args = parser.parse_args()
    create_test_h5ad(args.output, args.nobs, args.nvar, args.nnz_percent, args.col_shift, args.seed)


def create_test_h5ad(outfile, nobs, nvar, nnz_percent=100, apply_col_shift=False, seed=None):
    random.seed(seed)
    np.random.seed(seed)
    x = create_X_array(nobs, nvar, nnz_percent, apply_col_shift)
    obsm = {"X_random": np.random.rand(nobs, 2).astype(np.float32)}
    adata = anndata.AnnData(x, obsm=obsm)
    adata.write(outfile)


def create_X_array(nobs, nvar, nnz_percent, apply_col_shift):
    if nnz_percent < 100:
        array = scipy.sparse.random(nobs, nvar, nnz_percent * 0.01, dtype=np.float32, format="csc")
    else:
        array = np.random.rand(nobs, nvar).astype(np.float32)

    if apply_col_shift:
        col_shift = np.random.rand((nvar))
        array += col_shift

    return array


if __name__ == "__main__":
    main()
