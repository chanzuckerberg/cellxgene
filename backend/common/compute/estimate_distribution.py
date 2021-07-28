import numba
import concurrent.futures
import numpy as np
from scipy import sparse
from backend.common.constants import XApproximateDistribution


@numba.njit(error_model="numpy", nogil=True)
def min_max(arr: np.ndarray):
    """Return (min, max) values for the ndarray."""

    # initialize to first finite value in array.  Normally,
    # this will exit on the first value.
    for i in range(arr.size):
        min_val = max_val = arr[i]
        if np.isfinite(min_val):
            break

    # now find min/max, unrolled by two
    odd = arr.size % 2
    unrolled_loop_limit = arr.size - 1 if odd else arr.size
    i = 0
    while i < unrolled_loop_limit:
        x = arr[i]
        y = arr[i + 1]

        # ignore non-finites
        x = x if np.isfinite(x) else min_val
        y = y if np.isfinite(y) else min_val

        if x > y:
            x, y = y, x
        min_val = min(x, min_val)
        max_val = max(y, max_val)
        i += 2

    # handle the tail if any
    if odd:
        x = arr[arr.size - 1]

        # ignore non-finites
        x = x if np.isfinite(x) else min_val

        min_val = min(x, min_val)
        max_val = max(x, max_val)

    return min_val, max_val


def estimate_approximate_distribution(X) -> XApproximateDistribution:
    """
    Estimate the distribution (normal, count) of the X matrix.

    Currently this is based upon the assumption that scRNA-seq data is
    exponentially distributed in its raw (count) form, and when logged,
    any (max-min) range in excess of 24 is implies tens of millions of
    observations of a single feature and so is extremely unlikely.
    """
    if X.dtype.kind not in ["i", "u", "f"]:
        raise TypeError(f"Unsupported matrix dtype: {X.dtype.name}")

    if X.size == 0:
        # default for empty array
        return XApproximateDistribution.NORMAL

    if sparse.isspmatrix_csc(X) or sparse.isspmatrix_csr(X):
        Xdata = X.data
    elif type(X) is np.ndarray:
        Xdata = X.reshape(
            X.size,
        )
    else:
        raise TypeError(f"Unsupported matrix format: {str(type(X))}")

    CHUNKSIZE = 1 << 24
    if Xdata.size > CHUNKSIZE:
        min_val = max_val = Xdata[0]
        with concurrent.futures.ThreadPoolExecutor() as tp:
            for (_min, _max) in tp.map(min_max, [Xdata[i : i + CHUNKSIZE] for i in range(0, Xdata.size, CHUNKSIZE)]):
                min_val = min(_min, min_val)
                max_val = max(_max, max_val)

    else:
        min_val, max_val = min_max(Xdata)

    excess_range = (max_val - min_val) > 24
    return XApproximateDistribution.COUNT if excess_range else XApproximateDistribution.NORMAL
