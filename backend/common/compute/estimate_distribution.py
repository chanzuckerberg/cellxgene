import numba
import concurrent.futures
import numpy as np
from scipy import sparse
from backend.common.constants import XApproxDistribution


def dtype_limits(arr: np.ndarray):
    if arr.dtype.kind == "f":
        finf = np.finfo(arr.dtype)
        return (finf.min, finf.max)
    if arr.dtype.kind in ["u", "i"]:
        iinf = np.iinfo(arr.dtype)
        return (iinf.min, iinf.max)
    # error - should never occur
    raise TypeError(f"Unsupported matrix dtype: {arr.dtype.name}")


@numba.njit(fastmath=True, error_model="numpy", nogil=True)
def min_max(arr: np.ndarray, min_limit, max_limit):
    """Return (min, max) values for the ndarray."""
    n = arr.size
    odd = n % 2
    if odd:
        n -= 1

    max_val, min_val = min_limit, max_limit

    i = 0
    while i < n:
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

    if odd:
        x = arr[n]

        # ignore non-finites
        x = x if np.isfinite(x) else min_val

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
    if X.dtype not in ["i", "u", "f"]:
        raise TypeError(f"Unsupported matrix dtype: {X.dtype.name}")

    if X.size == 0:
        # default - empty array
        return XApproxDistribution.NORMAL

    if sparse.isspmatrix_csc(X) or sparse.isspmatrix_csr(X):
        Xdata = X.data
    elif type(X) is np.ndarray:
        Xdata = X.reshape(
            X.size,
        )
    else:
        raise TypeError(f"Unsupported matrix format: {str(type(X))}")

    min_limit, max_limit = dtype_limits(Xdata)

    CHUNKSIZE = 1 << 24
    if Xdata.size > CHUNKSIZE:
        min_val = max_val = Xdata[0]
        with concurrent.futures.ThreadPoolExecutor() as tp:
            for (_min, _max) in tp.map(
                lambda p: min_max(*p),
                [[Xdata[i : i + CHUNKSIZE], min_limit, max_limit] for i in range(0, Xdata.size, CHUNKSIZE)],
            ):
                min_val = min(_min, min_val)
                max_val = max(_max, max_val)

    else:
        min_val, max_val = min_max(Xdata, min_limit, max_limit)

    excess_range = (max_val - min_val) > 24
    return XApproxDistribution.COUNT if excess_range else XApproxDistribution.NORMAL
