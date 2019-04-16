# cellxgene roadmap

We are very exited for _cellxgene_ to become a valuable tool in collaborations
between computational biologists and experimental biologists working on
single-cell transcriptomics data. _cellxgene_ is in active development, and we
would love to include the community as we plan new features to work on. If you
have questions of feedback about this roadmap, please submit an issue on
GitHub.

Please note: this roadmap may change at any time.

*Last updated: April 11, 2019*

## what we are building now

In the near term, our goal is to enable teams of computational and experimental
biologists to explore and annotate their single-cell RNA-seq data.

There are 4 key features we plan to implement in the near term.

- Click install and launch
- Manual annotation workflows
- Toggle embeddings
- Gene information

### simple install and launch

The command line interface for installing and launching cellxgene is a barrier
for users who are not used to Python or using the command line. We plan to
support installation and launch of cellxgene on Mac and Windows. See
https://github.com/chanzuckerberg/cellxgene/issues/687 for more details.

### manual annotation workflows

The exploratory visualization that cellxgene offers is critical for manual
annotation workflows, especially in collaborative environments. We plan to
support manually annotate cells with labels (i.e., cell type or QC flags) for
downstream analysis. See https://github.com/chanzuckerberg/cellxgene/issues/524
for more details.

### toggle embeddings

While a single dataset may have multiple embeddings calculated (tSNE, umap, in
situ coordinates, trajectories, etc), cellxgene currently requires the user to select the
embedding to use in the main layout at launch. We plan to support letting users
toggle between any embedding present in a file from the cellxgene interface.
See https://github.com/chanzuckerberg/cellxgene/issues/594 for details.

### gene information

Differential expression returns only the names of genes, but no additional information
about gene metadata, function, or known associations. We plan to help users learn
more about genes they discover by exposing additional gene metadata. See
https://github.com/chanzuckerberg/cellxgene/issues/96 for details.
