"""Measure sparse differential expression on deterministic synthetic inputs.

Run in a fresh Python process from the checkout root, setting PYTHONPATH to ".".
Pass --module-source a saved copy of the original diffexp_generic.py to compare.
"""

import argparse
import importlib.util
import resource
import sys
import time

import numpy as np
from scipy import sparse

from server.common.compute import diffexp_generic
from server.common.constants import XApproximateDistribution


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--module-source")
    args = parser.parse_args()
    algorithm = diffexp_generic
    if args.module_source:
        spec = importlib.util.spec_from_file_location("benchmark_diffexp", args.module_source)
        algorithm = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(algorithm)

    rows, columns, entries = 64_000, 4096, 96
    rng = np.random.default_rng(87)
    ptr = np.arange(rows + 1, dtype=np.int32) * entries
    indices = rng.integers(columns, size=rows * entries, dtype=np.int32)
    values = rng.random(rows * entries, dtype=np.float32)
    matrix = sparse.csr_matrix((values, indices, ptr), shape=(rows, columns))
    matrix.sum_duplicates()
    a = np.arange(rows) < rows // 2

    class Adaptor:
        def get_X_array(self, obs_mask, _var_mask):
            return matrix[obs_mask]

        def get_X_approximate_distribution(self):
            return XApproximateDistribution.NORMAL

    # Time the calculation after warming the reducer for each implementation.
    algorithm.mean_var_n(sparse.csr_matrix([[0.0], [1.0]]))
    start = time.perf_counter()
    results = algorithm.diffexp_ttest(Adaptor(), a, ~a, top_n=25)
    elapsed = time.perf_counter() - start
    assert len(results["positive"]) == len(results["negative"]) == 25
    mib = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024**2 if sys.platform == "darwin" else 1024)
    print(f"{matrix.shape}, nnz={matrix.nnz}: {elapsed:.3f} s, peak RSS {mib:.1f} MiB")
    print("Top positive gene indices:", [row[0] for row in results["positive"]])


if __name__ == "__main__":
    main()
