import numba
import concurrent.futures
import numpy as np
from scipy import sparse
from backend.common.constants import XApproxDistribution


@numba.njit(fastmath=True, error_model="numpy", nogil=True)
def min_max(arr):
    """Return (min, max) values for the ndarray."""
    n = arr.size
    odd = n % 2
    if not odd:
        n -= 1
    max_val = min_val = arr[0]
    i = 1
    while i < n:
        x = arr[i]
        y = arr[i + 1]
        if x > y:
            x, y = y, x
        min_val = min(x, min_val)
        max_val = max(y, max_val)
        i += 2
    if not odd:
        x = arr[n]
        min_val = min(x, min_val)
        max_val = max(x, max_val)
    return min_val, max_val


def estimate_approximate_distribution(X) -> XApproxDistribution:
    """
    Estimate the distribution (normal, count) of the X matrix.

    Currently this is based upon the assumption that scRNA-seq data is
    exponentially distributed in its raw (count) form, and when logged,
    any (max-min) range in excess of 24 is implies tens of millions of
    observations of a single feature and so is extremely unlikely.
    """
    if sparse.isspmatrix_csc(X) or sparse.isspmatrix_csr(X):
        Xdata = X.data
    elif type(X) is np.ndarray:
        Xdata = X.reshape(
            X.size,
        )
    else:
        raise TypeError(f"Unsupported matrix type: {str(type(X))}")

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
    return XApproxDistribution.COUNT if excess_range else XApproxDistribution.NORMAL
