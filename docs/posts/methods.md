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
**N.B.: this implementation assumes normally distributed values on a linear scale.**

Currently, we use a [Welch's _t_-test](https://en.wikipedia.org/wiki/Welch%27s_t-test), which assumes that the two populations are normally distributed (but may have unequal variance). We use a two-sided test against the null hypothesis that the two populations have **equal** means (i.e., based on the magnitude of the difference in means, regardless of directionality). We also use the same variance overestimation correction as in `scanpy`.

We sort genes by `|t values|`, and return the top 15 that have a `log fold change >= 0.01`. Both the number of genes returned and the log fold change threshold can be [configured in the CLI](launch).
