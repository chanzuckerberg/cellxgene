import os
import concurrent.futures
from itertools import repeat
import numpy as np
from server.compute.diffexp_generic import diffexp_ttest_from_mean_var
from server.data_cxg.cxg_util import pack_selector

"""
See the comments in diffexp_generic for a description of this algorithm

This implementation runs directly in-process.  It is multi- threaded, but not particularly scalable.
Longer term, will likely move to a distributed framework for this.

There are currently no global throttles on simultaneous workers.
"""

# number of simultaneous workers, per HTTP request
MAX_WORKERS = 16


def diffexp_ttest(adaptor, maskA, maskB, top_n=8, diffexp_lfc_cutoff=0.01):
    row_selector_A = np.where(maskA)[0]
    row_selector_B = np.where(maskB)[0]
    nA = len(row_selector_A)
    nB = len(row_selector_B)
    matrix = adaptor.open_array("X")
    # mean, variance, N - calculate for both selections
    with MyThreadPoolExecutor(max_workers=2) as executor:
        A = executor.submit(mean_var, matrix, row_selector_A)
        B = executor.submit(mean_var, matrix, row_selector_B)
        meanA, varA = A.result()
        meanB, varB = B.result()

    return diffexp_ttest_from_mean_var(meanA, varA, nA, meanB, varB, nB, top_n, diffexp_lfc_cutoff)


class DispatchMixins:
    def parallel_dispatch(self, fn, *iterables):
        """
        Dispatch jobs via concurrent.futures.Executor.submit() and wait for their
        completion, return result via future.result().  Primary purpose is to throttle
        dispatch rate so that only 'MAX_JOBS_QUEUE_LENGTH' jobs are running at any
        given time, reducing overall memory footprint.
        """
        MAX_JOBS_QUEUE_LENGTH = self._max_workers + 4

        def submit_more_jobs(jobs, active_jobs):
            for job in jobs:
                future = self.submit(fn, *job)
                active_jobs[future] = job
                if len(active_jobs) >= MAX_JOBS_QUEUE_LENGTH:
                    return True
            return False

        def result_iterator(jobs):
            active_jobs = {}  # map of future -> args
            try:
                while submit_more_jobs(jobs, active_jobs) or len(active_jobs) > 0:
                    for future in concurrent.futures.as_completed(active_jobs.keys()):
                        job = active_jobs[future]
                        result = future.result()
                        # be careful to not retain dangling references
                        del active_jobs[future], future
                        yield (result, job)
            except Exception as e:
                print(str(e))
                raise
            finally:
                for future in active_jobs.keys():
                    future.cancel()

        return result_iterator(zip(*iterables))


class MyThreadPoolExecutor(concurrent.futures.ThreadPoolExecutor, DispatchMixins):
    pass


def _mean_var(matrix, row_selector, col_range):
    X = matrix.multi_index[row_selector, col_range[0] : col_range[1] - 1][""]
    n = X.shape[0]
    mean = np.sum(X, axis=0, dtype=np.float64) / n
    dfm = X - mean
    sumsq = np.sum(np.multiply(dfm, dfm, dtype=np.float64), axis=0)
    variance = sumsq / (n - 1)
    return (mean, variance)


def mean_var(matrix, row_selector):
    """
    row_selector: list of row indices
    """
    dtype = matrix.dtype
    rows, cols = matrix.shape
    tile_extent = [dim.tile for dim in matrix.schema.domain]

    dispatch_func = _mean_var
    row_selector = pack_selector(row_selector)

    # because all IO is done per-tile, and we are always dense and col-major,
    # use the tile column size as the partition.  Revisit partitioning if we
    # change the X layout, or start using a non-local execution environment
    # which may have other constraints.
    cols_per_partition = tile_extent[1]
    col_partitions = [(c, min(c + cols_per_partition, cols)) for c in range(0, cols, cols_per_partition)]

    max_workers = max(1, min(MAX_WORKERS, os.cpu_count()))  # throttle max_workers

    mean = np.zeros((cols,), dtype=np.float64)
    var = np.zeros((cols,), dtype=np.float64)
    dispatch_args = [repeat(matrix), repeat(row_selector), col_partitions]
    with MyThreadPoolExecutor(max_workers=max_workers) as exec:
        for result in exec.parallel_dispatch(dispatch_func, *dispatch_args):
            # returns tuple: (return_val, dispatch_args)
            m, v = result[0]
            cols = result[1][2]
            mean[cols[0] : cols[1]] += m
            var[cols[0] : cols[1]] += v
            del result, m, v

    return (mean.astype(dtype), var.astype(dtype))
