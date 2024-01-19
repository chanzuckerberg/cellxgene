import concurrent.futures

import numpy as np
from numba import jit
from scipy import stats

from server.common.constants import XApproximateDistribution
from server.common.errors import ComputeError


def diffexp_ttest(adaptor, setA, setB, top_n=8, diffexp_lfc_cutoff=0.01, selector_lists=False):
    """
    Return differential expression statistics for top N variables.

    Algorithm:
    - compute fold change
    - compute Welch's t-test statistic and pvalue (w/ Bonferroni correction)
    - return top N abs(logfoldchange) where lfc > diffexp_lfc_cutoff

    If there are not N which meet criteria, augment by removing the logfoldchange
    threshold requirement.

    Notes on alogrithm:
    - Welch's ttest provides basic statistics test.
      https://en.wikipedia.org/wiki/Welch%27s_t-test
    - p-values adjusted with Bonferroni correction.
      https://en.wikipedia.org/wiki/Bonferroni_correction

    :param adaptor: DataAdaptor instance
    :param setA: observation selector (mask or list) for set 1
    :param setB: observation selector (mask or list) for set 2
    :param top_n: number of variables to return stats for
    :param diffexp_lfc_cutoff: minimum
    :param selector_lists: if True, the selectors are presumed to be masks of length n_obs; else, lists of obs indices
    absolute value returning [ varindex, logfoldchange, pval, pval_adj ] for top N genes
    :return:  for top N genes, {"positive": for top N genes, [ varindex, foldchange, pval, pval_adj ],
              "negative": for top N genes, [ varindex, foldchange, pval, pval_adj ]}
    """
    matrix = adaptor.open_X_array()
    dtype = matrix.attr(0).dtype
    n_obs, cols = adaptor.get_shape()
    is_sparse = matrix.schema.sparse

    if selector_lists:
        assert 0 <= setA[0] < n_obs
        assert 0 <= setA[-1] < n_obs
        assert 0 <= setB[0] < n_obs
        assert 0 <= setB[-1] < n_obs
        row_selector_A = setA
        row_selector_B = setB
    else:
        assert len(setA) == len(setB) == n_obs
        row_selector_A = setA.nonzero()[0]
        row_selector_B = setB.nonzero()[0]

    mean_var_cnt_fn = mean_var_cnt_sparse if is_sparse else mean_var_cnt_dense
    with concurrent.futures.ThreadPoolExecutor() as tp:
        A = tp.submit(mean_var_cnt_fn, matrix, cols, row_selector_A)
        B = tp.submit(mean_var_cnt_fn, matrix, cols, row_selector_B)
        try:
            meanA, varA, nA = A.result()
            meanB, varB, nB = B.result()
        except Exception as e:
            raise ComputeError(str(e)) from None

    return diffexp_ttest_from_mean_var(
        meanA=meanA.astype(dtype),
        varA=varA.astype(dtype),
        nA=nA,
        meanB=meanB.astype(dtype),
        varB=varB.astype(dtype),
        nB=nB,
        top_n=top_n,
        diffexp_lfc_cutoff=diffexp_lfc_cutoff,
    )


def diffexp_ttest_from_mean_var(meanA, varA, nA, meanB, varB, nB, top_n, diffexp_lfc_cutoff):
    # IMPORTANT NOTE: this code assumes the data is normally distributed and/or already logged.
    n_var = meanA.shape[0]
    top_n = min(top_n, n_var)

    # variance / N
    vnA = varA / min(nA, nB)  # overestimate variance, would normally be nA
    vnB = varB / min(nA, nB)  # overestimate variance, would normally be nB
    sum_vn = vnA + vnB

    # degrees of freedom for Welch's t-test
    with np.errstate(divide="ignore", invalid="ignore"):
        dof = sum_vn**2 / (vnA**2 / (nA - 1) + vnB**2 / (nB - 1))
    dof[np.isnan(dof)] = 1

    # Welch's t-test score calculation
    with np.errstate(divide="ignore", invalid="ignore"):
        tscores = (meanA - meanB) / np.sqrt(sum_vn)
    tscores[np.isnan(tscores)] = 0

    # p-value
    pvals = stats.t.sf(np.abs(tscores), dof) * 2
    pvals_adj = pvals * n_var
    pvals_adj[pvals_adj > 1] = 1  # cap adjusted p-value at 1

    # log fold change. The data is normally distributed/logged, so just subtract the means.
    logfoldchanges = meanA - meanB

    stats_to_sort = tscores
    # find all with lfc > cutoff
    lfc_above_cutoff_idx = np.nonzero(np.abs(logfoldchanges) > diffexp_lfc_cutoff)[0]

    # derive sort order
    if lfc_above_cutoff_idx.shape[0] > top_n * 2:
        # partition top N
        rel_t_partition = np.argpartition(stats_to_sort[lfc_above_cutoff_idx], (top_n, -top_n))
        rel_t_partition_top_n = np.concatenate((rel_t_partition[-top_n:], rel_t_partition[:top_n]))
        t_partition = lfc_above_cutoff_idx[rel_t_partition_top_n]
        # sort the top N partition
        rel_sort_order = np.argsort(stats_to_sort[t_partition])[::-1]
        sort_order = t_partition[rel_sort_order]
    else:
        # partition and sort top N, ignoring lfc cutoff
        partition = np.argpartition(stats_to_sort, (top_n, -top_n))
        partition_top_n = np.concatenate((partition[-top_n:], partition[:top_n]))

        rel_sort_order = np.argsort(stats_to_sort[partition_top_n])[::-1]
        indices = np.indices(stats_to_sort.shape)[0]
        sort_order = indices[partition_top_n][rel_sort_order]

    # top n slice based upon sort order
    logfoldchanges_top_n = logfoldchanges[sort_order]
    pvals_top_n = pvals[sort_order]
    pvals_adj_top_n = pvals_adj[sort_order]

    # varIndex, logfoldchange, pval, pval_adj
    result = {
        "positive": [
            [sort_order[i], logfoldchanges_top_n[i], pvals_top_n[i], pvals_adj_top_n[i]] for i in range(top_n)
        ],
        "negative": [
            [sort_order[i], logfoldchanges_top_n[i], pvals_top_n[i], pvals_adj_top_n[i]]
            for i in range(-1, -1 - top_n, -1)
        ],
    }

    return result


def mean_var_n(X, X_approximate_distribution=XApproximateDistribution.NORMAL):
    """
    Two-pass variance calculation.  Numerically (more) stable
    than naive methods (and same method used by numpy.var())
    https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Two-pass
    """
    # fp_err_occurred is a flag indicating that a floating point error
    # occurred somewhere in our compute.  Used to trigger non-finite
    # number handling.
    fp_err_occurred = False

    def fp_err_set(err, flag):
        nonlocal fp_err_occurred
        fp_err_occurred = True

    with np.errstate(divide="call", invalid="call", call=fp_err_set):
        n = X.shape[0]
        # TODO: fix to pull distribution from cxg schema? Or is this redundant
        if X_approximate_distribution == XApproximateDistribution.COUNT:
            X = np.log1p(X)
        mean = X.mean(axis=0)
        dfm = X - mean
        sumsq = np.sum(np.multiply(dfm, dfm), axis=0)
        v = sumsq / (n - 1)

    if fp_err_occurred:
        mean[np.isfinite(mean) == False] = 0  # noqa: E712
        v[np.isfinite(v) == False] = 0  # noqa: E712
    else:
        mean[np.isnan(mean)] = 0
        v[np.isnan(v)] = 0

    return mean, v, n


def mean_var_cnt_dense(matrix, _, rows):
    xslc = matrix.multi_index[rows]
    return mean_var_n(xslc[""])


@jit(nopython=True, nogil=True, fastmath=True)
def _mean_var_sparse_accumulate(col_arr, val_arr, n, u, M2):
    """
    Incrementally accumulate mean and sum of square of distance from mean using
    Welford's online method.
    """
    for col, val in zip(col_arr, val_arr):
        u_prev = u[col]
        M2_prev = M2[col]
        n[col] += 1
        u[col] = u_prev + (val - u_prev) / n[col]
        M2[col] = M2_prev + (val - u_prev) * (val - u[col])


@jit(nopython=True, nogil=True, fastmath=True)
def _mean_var_sparse_finalize(n_rows, n_a, u_a, M2_a):
    """
    Finalize incremental values, acconting for missing elements (due to sparse input).
    Non-sparse and sparse combined using Chan's parallel adaptation of Welford's.
    The code assumes the sparse elements are all zero and ignores those terms.
    """
    n_b = n_rows - n_a
    delta = -u_a  # assumes u_b == 0
    u = (n_a * u_a) / n_rows
    M2 = M2_a + delta**2 * n_a * n_b / n_rows  # assumes M2_b == 0
    return u, M2


def mean_var_cnt_sparse(matrix, n_var, rows):
    query_iterator = matrix.query(order="U", return_incomplete=True).multi_index[rows]

    # accumulators, by gene (var) for n, u (mean) and M (sum of squares of difference from mean)
    n_rows = len(rows)
    n_a = np.zeros((n_var,), dtype=np.uint32)
    u_a = np.zeros((n_var,), dtype=np.float64)
    M2_a = np.zeros((n_var,), dtype=np.float64)
    for slc in query_iterator:
        _mean_var_sparse_accumulate(slc["var"], slc[""], n_a, u_a, M2_a)

    u, M2 = _mean_var_sparse_finalize(n_rows, n_a, u_a, M2_a)

    # compute variance
    var = M2 / max(1, (n_rows - 1))

    return u, var, n_rows
