---
layout: default
title: prepare
description: Preparing your data
---
# Data format requirements

If your data is in `h5ad` file (from the [`anndata`](https://anndata.readthedocs.io/en/latest/index.html) library) and meets the following requirements, you can go straight to `cellxgene launch`:

- Expression values (raw or normalized) in `anndata.X`
- At least one embedding (e.g., tSNE, UMAP) in `anndata.obsm`, specified with the prefix `X_` (e.g., by default scanpy stores UMAP coordinates in `anndata.obsm['X_umap']`)
- A unique identifier is required for each cell, which by default will be pulled from the `obs` DataFrame index. If the index is not unique or does not contain the cell ID, an alternative column can be specified with `--obs-names`
- A unique identifier is required for each gene, which by default will be pulled from the `var` DataFrame index. If the index is not unique or does not contain the gene ID, an alternative column can be specified with `--var-names`

#### What about R objects from seurat / bioconductor!?
We hear you! We'd also love to be able to ingest these files directly. This isn't currently possible, but in the meantime, you can use [sceasy](https://bioconda.github.io/recipes/r-sceasy/README.html) ([docs](https://cellgeni.readthedocs.io/en/latest/visualisations.html)) to convert to `h5ad`. Seurat also has some [handy conversion tools](https://satijalab.org/seurat/v3.0/conversion_vignette.html) that you can try out.

#### Can I use data hosted on the web somewhere?
Yes! You can launch from a URL instead of a filepath. The same data format requirements apply. Please see [here](launch) for more details.

# Data format options

#### Category colors
`cellxgene` will display [scanpy-style color
information](https://github.com/chanzuckerberg/cellxgene/issues/1152#issue-564361541)
for category-label pairs. An example of this format is shown below:

```
>>> category = "louvain"
>>> # colors stored in adata.uns must be matplotlib-compatible color information
>>> adata.uns[f"{category}_colors"]
array(['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#bcbd22'], dtype='<U7')
>>> # there must be a matching category in adata.obs
>>> category in adata.obs
True
```

To test that you've done this properly, check that for your given `category` the number of colors match the number of category values and that the second command below results in a mapping from categories to colors.

```
>>> len(adata.obs[category].cat.categories) == len(adata.uns[f"{category}_colors"])
True
>>> dict(zip(adata.obs[category].cat.categories, adata.uns[f"{category}_colors"]))
{'CD4 T cells': '#1f77b4', 'CD14+ Monocytes': '#ff7f0e', 'B cells': '#2ca02c', 'CD8 T cells': '#d62728', 'NK cells': '#9467bd', 'FCGR3A+ Monocytes': '#8c564b', 'Dendritic cells': '#e377c2', 'Megakaryocytes': '#bcbd22'}
```

You can disable this feature using the `--disable-custom-colors` flag for `cellxgene launch`. cellxgene will then chose colors from its standard color palettes.

# Using `cellxgene prepare`

If your data is in a different format, and/or you still need to perform dimensionality reduction and/or clustering, `cellxgene` can do that for you with the `prepare` command.

## What is `cellxgene prepare`?

`cellxgene prepare` offers an easy command line interface (CLI) to preliminarily wrangle your data into the required format for previewing it with cellxgene. It runs `scanpy` under the hood and can read in any format that is currently supported by `scanpy` (including mtx, loom, and more listed [in the scanpy documentation](https://scanpy.readthedocs.io/en/latest/api/index.html#reading)).

`prepare` uses scanpy to:

- Handle simple data normalization (from a [recipe](https://www.pydoc.io/pypi/scanpy-0.2.3/autoapi/preprocessing/recipes/index.html))
- Do basic preprocessing to run PCA and compute the neighbor graph
- Reduce dimensionality to generate embeddings
- Infer clusters

You can control which steps to run and their methods (when applicable), via the CLI. The CLI also includes options for computing QC metrics, enforcing matrix sparcity, specifying index names, and plotting output.

## What is cellxgene `prepare` _not_?

`cellxgene prepare` is not meant as a way to formally process or analyze your data. It's simply a utility for quickly wrangling your data into cellxgene-compatible format and computing a "vanilla" embedding so you can try out `cellxgene` and get a general sense of a dataset.

## Quickstart for `cellxgene prepare`
To add `cellxgene prepare` to your [cellxgene installation](install), run
`pip install cellxgene[prepare]`

Then run `prepare` on your data with:
```
cellxgene prepare dataset.h5ad --output=dataset-processed.h5ad
```

This will load the input data, perform PCA and nearest neighbor calculations, compute `UMAP` and `tSNE` embeddings and `louvain` cluster assignments, and save the results in a new file called `dataset-processed.h5ad` that can be loaded using `cellxgene launch`.

## Example usage

As a quick example, let's construct a command to use `prepare` to take a raw expression matrix and generate a processed `h5ad` ready to visualize with cellxgene.

We'll start off using the raw data from the pbmc3k dataset. This dataset is described [here](https://icb-scanpy.readthedocs-hosted.com/en/stable/api/scanpy.datasets.pbmc3k.html), and is available as part of the scanpy package. For this example, we'll assume this raw data is stored in a file called `pbmc3k-raw.h5ad`.

Our `prepare` command looks like this:

```
cellxgene prepare pbmc3k-raw.h5ad \
	--run-qc \                                  # (A)
	--recipe seurat \                           # (B)
	--layout tsne --layout umap \               # (C)
	--output pbmc3k-prepared.h5ad               # (D)
```

Let's look at what `prepare` is doing to our data, and how each step relates to the command above. You can see a walkthrough of what's going on under the hood for this example in [this notebook](https://github.com/chanzuckerberg/cellxgene-vignettes/blob/master/dataset-processing/pbmc3k-prepare-example.ipynb).

**(A) - Compute quality control metrics and store this in our `AnnData` object for later inspection**
**(B) - Normalize the expression matrix using a basic preprocessing recipe**
**(auto) - Do some preprocessing to run PCA and compute the neighbor graph**
**(auto) - Infer clusters with the Louvain algorithm and store these labels to visualize later**
**(C) - Compute and store UMAP and tSNE embeddings**
**(D) - Write results to file**

## Options for cellxgene `prepare`

**For the most up-to-date and comprehensive list of options, run `cellxgene prepare --help`**

`--embedding` controls which dimensionality reduction algorithm is applies to your data.
Options are `umap` and/or `tsne`. Defaults to both.

`--recipe` controls which normalization steps to apply to your data, based on one of the preprocessing `recipes` included with `scanpy`.
These recipes include steps like cell filtering and gene selection; see the `scanpy` [documentation](https://scanpy.readthedocs.io/en/latest/api/index.html#recipes) for more details.
Options are `none`, `seurat`, or `zheng17`. Defaults to `none`.

`--sparse` is a flag determines whether to enforce a sparse matrix. For large datasets, `prepare` can take a long time to run (a few minutes for datasets with 10-100k cells, up to an hour or more for datasets with >100k cells). If you want `prepare` to run faster we recommend using the `sparse` option.
If this flag is not included, default is `False`

`--skip-qc` by default, `cellxgene prepare` will compute quality control metrics (saved to `anndata.obs` and `anndata.var`) as described in the `scanpy` [documentation](https://scanpy.readthedocs.io/en/stable/api/scanpy.pp.calculate_qc_metrics.html). Pass this flag if you would like to skip this step.

`--make-obs-names-unique` / `--make-var-names-unique` determine whether to rename `obs` (cell) / `var` (gene) names, respectively, to be unique.
Default is `True`.

`--set-obs-names` controls which field in `anndata.obs` (cell metadata) is used as the _index_ for cells (e.g., a cell ID column).
Default is `anndata.obs.names`

`--set-var-names` controls which field in `anndata.var` (gene metadata) is used as the _index_ for genes.
Default is `anndata.var.names`

`--output` and `--overwrite` control where the processed data is saved.
