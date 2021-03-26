"""
Script to create a sparse dataset in CXG format based on an input dataset in CXG format.
The input dataset is not modified.
"""
import argparse
import os
import shutil
import sys

import tiledb

from backend.czi_hosted.common.utils.cxg_generation_utils import convert_ndarray_to_cxg_dense_array, \
    convert_matrix_to_cxg_array
from backend.czi_hosted.common.utils.matrix_utils import is_matrix_sparse, get_column_shift_encode_for_matrix


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="input cxg directory")
    parser.add_argument("output", help="output cxg directory")
    parser.add_argument("--overwrite", action="store_true", help="replace output cxg directory")
    parser.add_argument("--verbose", "-v", action="count", default=0, help="verbose output")
    parser.add_argument(
        "--sparse-threshold",
        "-s",
        type=float,
        default=5.0,  # default is 5% non-zero values
        help="The X array will be sparse if the percent of non-zeros falls below this value",
    )
    args = parser.parse_args()

    if os.path.exists(args.output):
        print("output dir exists:", args.output)
        if args.overwrite:
            print("output dir removed:", args.output)
            shutil.rmtree(args.output)
        else:
            print("use the overwrite option to remove the output directory")
            sys.exit(1)

    if not os.path.isdir(args.input):
        print("input is not a directory", args.input)
        sys.exit(1)

    shutil.copytree(args.input, args.output, ignore=shutil.ignore_patterns("X", "X_col_shift"))

    ctx = tiledb.Ctx(
        {
            "sm.num_reader_threads": 32,
            "sm.num_writer_threads": 32,
            "sm.consolidation.buffer_size": 1 * 1024 * 1024 * 1024,
        }
    )

    with tiledb.DenseArray(os.path.join(args.input, "X"), mode="r", ctx=ctx) as X_in:
        x_matrix_data = X_in[:, :]
        matrix_container = args.output

        is_sparse = is_matrix_sparse(x_matrix_data, args.sparse_threshold)
        if not is_sparse:
            col_shift = get_column_shift_encode_for_matrix(x_matrix_data, args.sparse_threshold)
            is_sparse = col_shift is not None
        else:
            col_shift = None

        if col_shift is not None:
            x_col_shift_name = f"{args.output}/X_col_shift"
            convert_ndarray_to_cxg_dense_array(x_col_shift_name, col_shift, ctx)
            tiledb.consolidate(matrix_container, ctx=ctx)
        if is_sparse:
            convert_matrix_to_cxg_array(matrix_container, x_matrix_data, is_sparse, ctx, col_shift)
            tiledb.consolidate(matrix_container, ctx=ctx)

    if not is_sparse:
        print("The array is not sparse, cleaning up, abort.")
        shutil.rmtree(args.output)
        sys.exit(1)


if __name__ == "__main__":
    main()
