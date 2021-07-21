import numpy as np
from scipy import sparse, stats
from backend.common.constants import XApproxDistribution


def diffexp_ttest(adaptor, maskA, maskB, top_n=8, diffexp_lfc_cutoff=0.01):
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
    :param maskA: observation selection mask for set 1
    :param maskB: observation selection mask for set 2
    :param top_n: number of variables to return stats for
    :param diffexp_lfc_cutoff: minimum
    absolute value returning [ varindex, logfoldchange, pval, pval_adj ] for top N genes
    :return:  for top N genes, {"positive": for top N genes, [ varindex, foldchange, pval, pval_adj ], "negative": for top N genes, [ varindex, foldchange, pval, pval_adj ]}
    """

    X_approx_distribution = adaptor.get_X_approx_distribution()
    dataA = adaptor.get_X_array(maskA, None)
    dataB = adaptor.get_X_array(maskB, None)

    # mean, variance, N - calculate for both selections
    meanA, vA, nA = mean_var_n(dataA, X_approx_distribution)
    meanB, vB, nB = mean_var_n(dataB, X_approx_distribution)
    res = diffexp_ttest_from_mean_var(meanA, vA, nA, meanB, vB, nB, top_n, diffexp_lfc_cutoff)

    return res


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
        dof = sum_vn ** 2 / (vnA ** 2 / (nA - 1) + vnB ** 2 / (nB - 1))
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


# Convenience function which handles sparse data
def mean_var_n(X, X_approx_distribution=XApproxDistribution.NORMAL):
    """
    Two-pass variance calculation.  Numerically (more) stable
    than naive methods (and same method used by numpy.var())
    https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Two-pass
    """
    # fp_err_occurred is a flag indicating that a floating point error
    # occured somewhere in our compute.  Used to trigger non-finite
    # number handling.
    fp_err_occurred = False

    def fp_err_set(err, flag):
        nonlocal fp_err_occurred
        fp_err_occurred = True

    with np.errstate(divide="call", invalid="call", call=fp_err_set):
        n = X.shape[0]
        if sparse.issparse(X):
            if X_approx_distribution == XApproxDistribution.EXPONENTIAL:
                X = X.log1p()
            mean = X.mean(axis=0).A1
            dfm = X - mean
            sumsq = np.sum(np.multiply(dfm, dfm), axis=0).A1
            v = sumsq / (n - 1)
        else:
            if X_approx_distribution == XApproxDistribution.EXPONENTIAL:
                X = np.log1p(X)
            mean = X.mean(axis=0)
            dfm = X - mean
            sumsq = np.sum(np.multiply(dfm, dfm), axis=0)
            v = sumsq / (n - 1)

    # AnnData does not guarantee that operations on a view of X will
    # return an ndarray, so force the cast if it wasn't done for us.
    if type(mean) is not np.ndarray:
        mean = mean.toarray()
    if type(v) is not np.ndarray:
        v = v.toarray()

    if fp_err_occurred:
        mean[np.isfinite(mean) == False] = 0  # noqa: E712
        v[np.isfinite(v) == False] = 0  # noqa: E712
    else:
        mean[np.isnan(mean)] = 0
        v[np.isnan(v)] = 0

    return mean, v, n
