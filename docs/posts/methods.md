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

Currently, we use a [Welch's _t_-test](https://en.wikipedia.org/wiki/Welch%27s_t-test) implementation, including the same variance overestimation correction as used in `scanpy`. We sort the `tscore` to identify the top 15 genes, and then filter to remove any that fall below a cutoff log fold change value, which can help remove spurious test results. The default threshold is `0.01` and can be changed using the option `--diffexp-lfc-cutoff`.
