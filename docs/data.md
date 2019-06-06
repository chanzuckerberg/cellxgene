---
layout: default
title: data
description: Data
---

# Data format requirements

# Data vignette: about cellxgene prepare

#### What is `cellxgene prepare`?

`prepare` offers an easy command line interface (CLI) to wrangle your data into the required format for use with `cellxgene`.

#### What is `cellxgene prepare` _not_?

`cellxgene prepare` is not meant as a way to formally process your data. It's simply a utility for quickly wrangling your data into cellxgene-compatible [format]() and computing a "vanilla" embedding so you can try out `cellxgene` and get a general sense of a dataset.

#### What input formats does it accept?

Currently, we accept `h5ad` and `loom` files, as well as `10x` directories, and are hoping to accept more formats in the future.  
While we'd like to support quick conversion from seurat and bioconductor, these packages don't currently output a python-parseable intermediate file type. In the meantime, you might check out the [converters](https://satijalab.org/seurat/v3.0/conversion_vignette.html) that are under early development.

#### What can `cellxgene prepare` do?

`prepare` uses scanpy to handle data normalization (from a [recipe](https://www.pydoc.io/pypi/scanpy-0.2.3/autoapi/preprocessing/recipes/index.html)); basic preprocessing to run PCA and compute the neighbor graph; infer clusters; and reduce dimensionality to generate embeddings. You can control which steps to run, and their methods and parameters, via the CLI. The CLI also includes options for computing QC metrics, enforcing matrix sparcity, specifying index names, and plotting output.

**To see a full list of available arguments and options, run `cellxgene prepare --help`.**

#### How do I use `cellxgene prepare`?

As a quick example, let's construct a command to use `prepare` to take a raw expression matrix and generate a processed `h5ad` ready to visualize with cellxgene. You can see a walkthrough of what `prepare` is doing under the hood for this example in [this notebook](????).

We'll start off using the raw data from the pbmc3k dataset. This dataset is described [here](https://icb-scanpy.readthedocs-hosted.com/en/stable/api/scanpy.datasets.pbmc3k.html), and is available as part of the scanpy API. For this example, we'll assume this raw data is stored in a file called `pbmc3k-raw.h5ad`.

We can start to compose our command like so:  
`cellxgene prepare pbmc3k-raw.h5ad`

We'll then add on some additional arguments to specify how we want `prepare` to handle our data:

**1 - Compute quality control metrics and store this in our `AnnData` object for later inspection**  
`cellxgene prepare pbmc3k-raw.h5ad --run-qc`

**2 - Normalize the expression matrix using a basic preprocessing recipe**  
`cellxgene prepare pbmc3k-raw.h5ad --run-qc --recipe seurat`

**3 - Do some preprocessing to run PCA and compute the neighbor graph**  
These are run automatically as they are required for the downstream steps, so we don't need an argument for this.

**4 - Infer clusters with the Louvain algorithm and store these labels to visualize later**  
`cellxgene prepare` always runs louvain clustering, so we don't need to pass an argument for this, either.

**5 - Compute and store umap and tsne embeddings**  
`cellxgene prepare --data pbmc3k-raw.h5ad --run-qc --recipe seurat --layout tsne --layout umap`

**6 - Write our results to file**  
`cellxgene prepare pbmc3k-raw.h5ad --run-qc --recipe seurat --layout tsne --layout umap --output pbmc3k-prepared.h5ad`

# Example datasets to use with cellxgene

Each dataset was preprocessed for use with cellxgene, following the original QC and preprocessing methods as closely as possible.  
All datasets include several embeddings and labels from multiple clustering methods (see linked notebooks).

**To download and use these datasets, run:**  
`curl -O [downloadURL]`  
`unzip [filename.zip]`  
`cellxgene launch [filename.h5ad] --open`

### Peripheral blood mononuclear cells

Healthy human PBMCs (10X).

- Source: [10X genomics](https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/pbmc3k)
- Cells: 2,638
- File size: 19MB
- [Raw data](http://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz)
- [Processing](https://github.com/chanzuckerberg/cellxgene-vignettes/blob/master/dataset-processing/pbmc3k-processing.ipynb)
- Download URL: ``

### Tabula muris

20 organs and tissues from healthy mice (Smart-Seq2).  
Rich metadata and annotations.

- Source: [CZBiohub](https://www.biorxiv.org/content/10.1101/237446v2)
- Cells: 45,423
- File size: 174MB
- [Raw data](https://figshare.com/projects/Tabula_Muris_Transcriptomic_characterization_of_20_organs_and_tissues_from_Mus_musculus_at_single_cell_resolution/27733)
- [Processing](https://github.com/chanzuckerberg/cellxgene-vignettes/blob/master/dataset-processing/tabula-muris-processing.ipynb)
- Download URL: ``

### Tabula muris senis

22 organs and tissues from healthy mice at ages 3mo, 18mo, 21mo, and 24mo (Smart-Seq2).  
Rich metadata and annotations.

- Source: [CZBiohub]()
- Cells: 81,478
- File size: 3.9GB
- [Raw data]()
- [Processing]()
- Download URL: ``
