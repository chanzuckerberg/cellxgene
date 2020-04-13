import concurrent.futures
import numpy as np
from server.compute.diffexp_generic import diffexp_ttest_from_mean_var, mean_var_n
from server.data_cxg.cxg_util import pack_selector_from_indices
from server.common.errors import ComputeError

"""
See the comments in diffexp_generic for a description of this algorithm

This implementation runs directly in-process.  It is multi- threaded, but not particularly scalable.
Longer term, will likely move to a distributed framework for this.

There are currently no global throttles on simultaneous workers.
"""

diffexp_thread_executor = None
max_workers = None
target_workunit = None


def set_config(config_max_workers, config_target_workunit):
    global max_workers
    global target_workunit
    max_workers = config_max_workers
    target_workunit = config_target_workunit


def get_thread_executor():
    global diffexp_thread_executor
    if diffexp_thread_executor is None:
        diffexp_thread_executor = concurrent.futures.ThreadPoolExecutor(max_workers=max_workers)
    return diffexp_thread_executor


def diffexp_ttest(adaptor, maskA, maskB, top_n=8, diffexp_lfc_cutoff=0.01):
    row_selector_A = np.where(maskA)[0]
    row_selector_B = np.where(maskB)[0]
    nA = len(row_selector_A)
    nB = len(row_selector_B)
    matrix = adaptor.open_array("X")

    cols = matrix.shape[1]
    tile_extent = [dim.tile for dim in matrix.schema.domain]

    # The rows from both row_selector_A and row_selector_B are gathered at the
    # same time, then the mean and variance are computed by subsetting on that
    # combined submatrix.  Combining the gather reduces number of requests/bandwidth
    # to the data source.
    row_selector_AB = np.union1d(row_selector_A, row_selector_B)
    row_selector_A_in_AB = np.in1d(row_selector_AB, row_selector_A, assume_unique=True)
    row_selector_B_in_AB = np.in1d(row_selector_AB, row_selector_B, assume_unique=True)
    row_selector_AB = pack_selector_from_indices(row_selector_AB)

    # because all IO is done per-tile, and we are always dense and col-major,
    # use the tile column size as the unit of partition.  Possibly access
    # more than one column tile at a time based on the target_workunit.
    # Revisit partitioning if we change the X layout, or start using a non-local execution environment
    # which may have other constraints.

    # TODO: If the number of row selections is large enough, then the cells_per_coltile will exceed
    # the target_workunit.  A potential improvement would be to partition by both columns and rows.
    # However partitioning the rows is slightly more complex due to the arbitrary distribution
    # of row selections that are passed into this algorithm.
    cells_per_coltile = (nA + nB) * tile_extent[1]
    cols_per_partition = max(1, int(target_workunit / cells_per_coltile)) * tile_extent[1]
    col_partitions = [(c, min(c + cols_per_partition, cols)) for c in range(0, cols, cols_per_partition)]

    meanA = np.zeros((cols,), dtype=np.float64)
    varA = np.zeros((cols,), dtype=np.float64)
    meanB = np.zeros((cols,), dtype=np.float64)
    varB = np.zeros((cols,), dtype=np.float64)

    executor = get_thread_executor()
    futures = []
    for cols in col_partitions:
        futures.append(
            executor.submit(_mean_var_ab, matrix, row_selector_AB, row_selector_A_in_AB, row_selector_B_in_AB, cols)
        )

    for future in futures:
        # returns tuple: (meanA, varA, meanB, varB, cols)
        try:
            result = future.result()
            part_meanA, part_varA, part_meanB, part_varB, cols = result
            meanA[cols[0]: cols[1]] += part_meanA
            varA[cols[0]: cols[1]] += part_varA
            meanB[cols[0]: cols[1]] += part_meanB
            varB[cols[0]: cols[1]] += part_varB
        except Exception as e:
            for future in futures:
                future.cancel()
            raise ComputeError(str(e))

    r = diffexp_ttest_from_mean_var(meanA, varA, nA, meanB, varB, nB, top_n, diffexp_lfc_cutoff)
    return r


def _mean_var_ab(matrix, row_selector_AB, row_selector_A_in_AB, row_selector_B_in_AB, col_range):
    X = matrix.multi_index[row_selector_AB, col_range[0] : col_range[1] - 1][""]
    meanA, varA, n = mean_var_n(X[row_selector_A_in_AB])
    meanB, varB, n = mean_var_n(X[row_selector_B_in_AB])
    return (meanA, varA, meanB, varB, col_range)
