---
layout: default
title: data
description: Data
---

# Using `cellxgene prepare`

#### What is `cellxgene prepare`?

`prepare` offers an easy command line interface (CLI) to preliminarily wrangle your data into the required format for previewing it with `cellxgene`.

#### What is `cellxgene prepare` _not_?

`cellxgene prepare` is not meant as a way to formally process or analyze your data. It's simply a utility for quickly wrangling your data into cellxgene-compatible format and computing a "vanilla" embedding so you can try out `cellxgene` and get a general sense of a dataset.

#### What input formats does it accept?

Currently, we accept `h5ad` and `loom` files, as well as `10x` directories, and are hoping to accept more formats in the future.

While we'd like to support quick conversion from seurat and bioconductor, these packages don't currently output a python-parseable intermediate file type. In the meantime, you might check out the [converters](https://satijalab.org/seurat/v3.0/conversion_vignette.html) that are under early development.

#### What can `cellxgene prepare` do?

`prepare` uses scanpy to:

- Handle simple data normalization (from a [recipe](https://www.pydoc.io/pypi/scanpy-0.2.3/autoapi/preprocessing/recipes/index.html))
- Do basic preprocessing to run PCA and compute the neighbor graph
- Infer clusters
- Reduce dimensionality to generate embeddings.  
  You can control which steps to run and their methods (when applicable), via the CLI. The CLI also includes options for computing QC metrics, enforcing matrix sparcity, specifying index names, and plotting output.

**To see a full list of available arguments and options, run `cellxgene prepare --help`.**

#### How do I use `cellxgene prepare`?

As a quick example, let's construct a command to use `prepare` to take a raw expression matrix and generate a processed `h5ad` ready to visualize with cellxgene.

We'll start off using the raw data from the pbmc3k dataset. This dataset is described [here](https://icb-scanpy.readthedocs-hosted.com/en/stable/api/scanpy.datasets.pbmc3k.html), and is available as part of the scanpy API. For this example, we'll assume this raw data is stored in a file called `pbmc3k-raw.h5ad`.

Our `prepare` compose our command looks like this:

```
cellxgene prepare pbmc3k-raw.h5ad \
	--run-qc \                                  # (A)
	-- --recipe seurat \                        # (B)
	--layout tsne --layout umap \               # (C)
	--output pbmc3k-prepared.h5ad               # (D)
```

Let's look at what `prepare` is doing to our data, and how each step relates to the command above. You can see a walkthrough of what's going on under the hood for this example in [this notebook](https://github.com/chanzuckerberg/cellxgene-vignettes/blob/master/dataset-processing/pbmc3k-prepare-example.ipynb).

**1 - Compute quality control metrics and store this in our `AnnData` object for later inspection (A)**  
**2 - Normalize the expression matrix using a basic preprocessing recipe (B)**  
**3 - Do some preprocessing to run PCA and compute the neighbor graph (auto)**  
**4 - Infer clusters with the Louvain algorithm and store these labels to visualize later (auto)**  
**5 - Compute and store umap and tsne embeddings (C)**  
**6 - Write results to file (D)**

# Example datasets to use with cellxgene

**To download and use these datasets, run:**
`curl -O [URL]`
`unzip [filename.zip]`
`cellxgene launch [filename.h5ad] --open`

### Peripheral blood mononuclear cells

Healthy human PBMCs (10X).

- Source: [10X genomics](https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/pbmc3k)
- Cells: 2,638
- File size: 19MB
- [Raw data](http://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz)
- [Processing](https://github.com/chanzuckerberg/cellxgene-vignettes/blob/master/dataset-processing/pbmc3k-processing.ipynb)
- Download: `curl -O https://cellxgene-example-data.czi.technology/pbmc3k.h5ad.zip`

### Tabula muris

20 organs and tissues from healthy mice (Smart-Seq2).  
Rich metadata and annotations.

- Source: [bioRxiv, CZBiohub](https://www.biorxiv.org/content/10.1101/237446v2)
- Cells: 45,423
- File size: 174MB
- [Raw data](https://figshare.com/projects/Tabula_Muris_Transcriptomic_characterization_of_20_organs_and_tissues_from_Mus_musculus_at_single_cell_resolution/27733)
- [Processing](https://github.com/chanzuckerberg/cellxgene-vignettes/blob/master/dataset-processing/tabula-muris-processing.ipynb)
- Download: `curl -O https://cellxgene-example-data.czi.technology/tabula-muris.h5ad.zip`

### Tabula muris senis

22 organs and tissues from healthy mice at ages 3mo, 18mo, 21mo, and 24mo (Smart-Seq2).  
Rich metadata and annotations.

- Source: [bioRxiv, CZBiohub](https://www.biorxiv.org/content/10.1101/661728v1)
- Cells: 81,478
- File size: 3.9GB
- Raw data [geo link coming soon!]
- [Processing](https://www.biorxiv.org/content/10.1101/661728v1)
- Download: `curl -O https://cellxgene-example-data.czi.technology/tabula-muris-senis.h5ad.zip`
