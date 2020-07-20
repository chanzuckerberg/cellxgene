# Methods

## Data structure: _anndata_ fields used for visualization

### Expression data  

Gene expression values are pulled from `anndata.X`. These feed into the histograms, scatterplot, colorscale, and differential expression calculations. We're [working on ways](https://github.com/chanzuckerberg/cellxgene/issues/689) to incorporate `anndata.raw` and other `anndata.layers`!

### Metadata  

Categorical (e.g., cluster labels) and continuous (e.g., pseudotime) metadata are pulled from `anndata.obs`. Any column added here will be available for visualization in cellxgene. You can also [create new categorical annotations](annotations) within the application.

### Embeddings

cellxgene looks for embeddings (e.g., tSNE, UMAP, PCA, spatial coordinates) in `anndata.obsm`. These fields must follow the scanpy convention of starting with `X_`, e.g., `anndata.obsm['X_umap']`. If an embedding has more than two components, the first two will be used for visualization.

## Differential expression

We're actively working on how to improve differential expression within the app.
**N.B.: the [current implementation](https://github.com/chanzuckerberg/cellxgene/blob/main/server/app/scanpy_engine/diffexp.py#L40) assumes normally distributed values on a linear scale.**

Currently, we use a [Welch's _t_-test](https://en.wikipedia.org/wiki/Welch%27s_t-test), which assumes that the two populations are each normally distributed, but may have unequal variance. We use a two-sided t-test against the null hypothesis that the two populations have **equal** means. P-values are adjusted with the [Bonferroni corrrection](https://en.wikipedia.org/wiki/Bonferroni_correction).

To help avoid spurious results, we use the log fold change to filter genes, retaining those where `|log2( mean(set1) / mean(set2) )| > 0.01`; this threshold can be configured with the [`--diffexp-lfc-cutoff`](launch) command. We then sort genes by their associated `|t value|` and return the top 15 genes.
