
import numpy as np
from scipy import sparse, stats


# Convenience function which handles sparse data
def _mean_var_n(X):
    """
    Two-pass variance calculation.  Numerically (more) stable
    than naive methods (and same method used by numpy.var())
    https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Two-pass
    """
    n = X.shape[0]
    if sparse.issparse(X):
        mean = X.mean(axis=0).A1
        dfm = X - mean
        sumsq = np.sum(np.multiply(dfm, dfm), axis=0).A1
        v = sumsq / (n - 1)
    else:
        mean = X.mean(axis=0)
        dfm = X - mean
        sumsq = np.sum(np.multiply(dfm, dfm), axis=0)
        v = sumsq / (n - 1)

    return mean, v, n


def diffexp_ttest(adata, maskA, maskB, top_n=8):
    """
    Return differential expression statistics for top N variables, sorted by
    t statistic.   Implemented as a unequal variance t-test.

    :param adata: anndata dataframe
    :param maskA: observation selection mask for set 1
    :param maskB: observation selection mask for set 2
    :param top_n: number of variables to return stats for
    :return:  for top N genes, [ varindex, logfoldchange, pval, pval_adj ]
    """

    # mean, variance, N
    meanA, vA, nA = _mean_var_n(adata._X[maskA])
    meanB, vB, nB = _mean_var_n(adata._X[maskB])

    # variance / N
    vnA = vA / nA
    vnB = vB / nB
    sum_vn = vnA + vnB

    # degrees of freedom for Welch's t-test
    with np.errstate(divide='ignore', invalid='ignore'):
        dof = sum_vn**2 / (vnA**2 / (nA - 1) + vnB**2 / (nB - 1))
    dof[np.isnan(dof)] = 1

    # Welch's t-test score calculation
    with np.errstate(divide='ignore', invalid='ignore'):
        tscores = (meanA - meanB) / np.sqrt(sum_vn)
    tscores[np.isnan(tscores)] = 0

    # p-value
    pvals = stats.t.sf(np.abs(tscores), dof) * 2
    pvals_adj = pvals * adata._X.shape[1]

    # logfoldchanges: log2(meanA / meanB)
    logfoldchanges = np.log2(np.abs((meanA + 1e-9) / (meanB + 1e-9)))

    # top n sort
    stats_to_sort = np.abs(tscores)
    partition = np.argpartition(stats_to_sort, -top_n)[-top_n:]
    rel_sort_order = np.argsort(stats_to_sort[partition])[::-1]
    vars_indices = np.arange(adata.n_vars, dtype=int)
    sort_order = vars_indices[partition][rel_sort_order]

    # top n slice
    logfoldchanges_top_n = logfoldchanges[sort_order]
    pvals_top_n = pvals[sort_order]
    pvals_adj_top_n = pvals_adj[sort_order]

    # varIndex, logfoldchange, pval, pval_adj
    result = [[sort_order[i],
               logfoldchanges_top_n[i],
               pvals_top_n[i],
               pvals_adj_top_n[i]] for i in range(top_n)]
    return result
