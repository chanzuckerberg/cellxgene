---
layout: default
title: roadmap
description: Roadmap
---

# Roadmap

cellxgene makes it easier for biologists to collaboratively explore and understand their single-cell RNA-seq data.
In the near term, we are focused on continuing to enable fast, interactive exploration of single-cell data, supporting collaborative workflows in single-cell analysis, and improving user support.
If you have questions or feedback about this roadmap, please submit an issue on GitHub.
Please note: this roadmap is subject to change.

Last updated: June 25, 2019

## Fast, interactive exploration of single-cell data

### Exposing Relationships Between Metadata and Data
Biologists need to understand how variables (stored in metadata) are associated with one another and how they relate to changes in gene expression.
Building upon visualization features that reveal categorical metadata relationships (cluster occupancy) and gene expression relationships (scatterplot), we plan to add exploratory visualization components that enable investigation of relationships between metadata and gene expression.
See [issue #616](https://github.com/chanzuckerberg/cellxgene/issues/616) for more details.

### Contextualizing Genes
While exploring a transcriptomics dataset, scientists need to understand the biological context of genes.
This context may be provided by user-defined gene metadata or publicly available gene databases.
We plan to support augmenting gene names with additional information that is useful to biologists.
See [issue #96](https://github.com/chanzuckerberg/cellxgene/issues/96) for more detail.

## Support collaborative workflows in single-cell analysis

### Manual Annotations
cellxgene offers exploratory visualizations that are critical for manual annotation workflows, especially in collaborative environments.
We plan to support manually annotating cells with labels (i.e., cell type or QC flags), and their easy export for downstream analysis.
See [issue #524](https://github.com/chanzuckerberg/cellxgene/issues/524) for more details.

### Simple Click to Launch [Paused]

Many biologists prefer not to interact with the command line and need an OS-native experience when using cellxgene.
We plan to implement a point-and-click installation and launch experience so that users can easily load data into cellxgene.
See [issue #687](https://github.com/chanzuckerberg/cellxgene/issues/687) for details.

### Python API [Paused]
For computational biologists, saving h5ad files then loading them into cellxgene is a point of friction.
We plan to support importing cellxgene as a Python package so that users can launch cellxgene directly from an interactive environment (such as Jupyter, IPython, or Spyder), and pass data to and from the cellxgene UI.

## Improving user support

### Improved documentation
cellxgene has some specific expectations about how data is stored.
We want to ensure that new users can get started easily and learn how to use cellxgene with their own data.
We plan to improve documentation on getting started, installation, data, and contributing.
See [issue #533](https://github.com/chanzuckerberg/cellxgene/issues/533) for more details.
