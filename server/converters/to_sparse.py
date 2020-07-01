"""
Script to create a sparse dataset in CXG format based on an input dataset in CXG format.
The input dataset is not modified.
"""
import os
import shutil
import tiledb
import argparse
import sys
import server.converters.cxgtool as cxgtool


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
        is_sparse = cxgtool.save_X(args.output, X_in, ctx, args.sparse_threshold, expect_sparse=True)

    if is_sparse is False:
        print("The array is not sparse, cleaning up, abort.")
        shutil.rmtree(args.output)
        sys.exit(1)


if __name__ == "__main__":
    main()
