---
layout: default
title: data
description: Data
---

# Data vignette: how to use cellxgene prepare

#### coming soon!

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
