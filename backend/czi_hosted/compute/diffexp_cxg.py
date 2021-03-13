import concurrent.futures
import numpy as np

from numba import jit

from backend.czi_hosted.data_cxg.cxg_util import pack_selector_from_indices
from backend.czi_hosted.compute.diffexp_generic import diffexp_ttest_from_mean_var, mean_var_n
from backend.czi_hosted.common.errors import ComputeError

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

    matrix = adaptor.open_array("X")
    row_selector_A = np.where(maskA)[0]
    row_selector_B = np.where(maskB)[0]
    nA = len(row_selector_A)
    nB = len(row_selector_B)

    dtype = matrix.dtype
    cols = matrix.shape[1]
    tile_extent = [dim.tile for dim in matrix.schema.domain]

    is_sparse = matrix.schema.sparse

    if is_sparse:
        row_selector_A = pack_selector_from_indices(row_selector_A)
        row_selector_B = pack_selector_from_indices(row_selector_B)
    else:
        # The rows from both row_selector_A and row_selector_B are gathered at the
        # same time, then the mean and variance are computed by subsetting on that
        # combined submatrix.  Combining the gather reduces number of requests/bandwidth
        # to the data source.
        row_selector_AB = np.union1d(row_selector_A, row_selector_B)
        row_selector_A_in_AB = np.in1d(row_selector_AB, row_selector_A, assume_unique=True)
        row_selector_B_in_AB = np.in1d(row_selector_AB, row_selector_B, assume_unique=True)
        row_selector_AB = pack_selector_from_indices(row_selector_AB)

    # because all IO is done per-tile, and we are always col-major,
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

    if is_sparse:
        for cols in col_partitions:
            futures.append(executor.submit(_mean_var_sparse_ab, matrix, row_selector_A, nA, row_selector_B, nB, cols))
    else:
        for cols in col_partitions:
            futures.append(
                executor.submit(_mean_var_ab, matrix, row_selector_AB, row_selector_A_in_AB, row_selector_B_in_AB, cols)
            )

    for future in futures:
        # returns tuple: (meanA, varA, meanB, varB, cols)
        try:
            result = future.result()
            part_meanA, part_varA, part_meanB, part_varB, cols = result
            meanA[cols[0] : cols[1]] += part_meanA
            varA[cols[0] : cols[1]] += part_varA
            meanB[cols[0] : cols[1]] += part_meanB
            varB[cols[0] : cols[1]] += part_varB
        except Exception as e:
            for future in futures:
                future.cancel()
            raise ComputeError(str(e))

    if is_sparse:
        if adaptor.has_array("X_col_shift"):
            X_col_shift = adaptor.open_array("X_col_shift")[:]
            meanA += X_col_shift
            meanB += X_col_shift

    r = diffexp_ttest_from_mean_var(
        meanA.astype(dtype),
        varA.astype(dtype),
        nA,
        meanB.astype(dtype),
        varB.astype(dtype),
        nB,
        top_n,
        diffexp_lfc_cutoff,
    )

    return r


def _mean_var_ab(matrix, row_selector_AB, row_selector_A_in_AB, row_selector_B_in_AB, col_range):
    X = matrix.multi_index[row_selector_AB, col_range[0] : col_range[1] - 1][""]
    meanA, varA, n = mean_var_n(X[row_selector_A_in_AB])
    meanB, varB, n = mean_var_n(X[row_selector_B_in_AB])
    return (meanA, varA, meanB, varB, col_range)


def _mean_var_sparse_ab(matrix, row_selector_A, nrows_A, row_selector_B, nrows_B, col_range):
    meanA, varA = _mean_var_sparse(matrix, row_selector_A, nrows_A, col_range)
    meanB, varB = _mean_var_sparse(matrix, row_selector_B, nrows_B, col_range)
    return (meanA, varA, meanB, varB, col_range)


@jit(nopython=True)
def _mean_var_sparse_numba(x, var, nrows, ncols):
    """Kernel to compute the mean and variance.  It was not clear if this function
    could be written using numpy, thus avoiding the loops.  Therefore numba is
    used here to speed things up.  With numba, this function takes a negligible amount
    of time compared to reading in the sparse matrix"""
    mean = np.zeros((ncols,), dtype=np.float64)
    for col, val in zip(var, x):
        mean[col] += val
    mean /= nrows

    # optimize the sumsq computation.
    # since most entries in a sparse matrix are 0, then start by assuming
    # all values are 0, so fill the sumsq array with nrows * (0 - mean)**2.
    # as non-zero values are encountered, subtract off the (mean*mean) value
    # and replace with (val-mean)**2.  Simplifying the expression
    # gives the following code.
    sumsq = nrows * np.multiply(mean, mean)
    for col, val in zip(var, x):
        sumsq[col] += val * (val - 2 * mean[col])
    v = sumsq / (nrows - 1)
    return mean, v


def _mean_var_sparse(matrix, selector, nrows, col_range):
    data = matrix.multi_index[selector, col_range[0] : col_range[1] - 1]
    x = data[""]

    # tiledb < 0.6.0 and >= 0.6.0 have slightly different interfaces.
    # the following takes care of both cases:
    #   older:  data["coords]["var"]
    #   newer:  data["var"]
    var = data.get("coords", data)["var"]

    # shift the column indices to start at 0, this
    # will become the index into the mean and var arrays.
    var -= col_range[0]

    fp_err_occurred = False

    def fp_err_set(err, flag):
        nonlocal fp_err_occurred
        fp_err_occurred = True

    ncols = col_range[1] - col_range[0]
    with np.errstate(divide="call", invalid="call", call=fp_err_set):
        mean, v = _mean_var_sparse_numba(x, var, nrows, ncols)

    if fp_err_occurred:
        mean[np.isfinite(mean) == False] = 0  # noqa: E712
        v[np.isfinite(v) == False] = 0  # noqa: E712
    else:
        mean[np.isnan(mean)] = 0
        v[np.isnan(v)] = 0

    return mean, v
