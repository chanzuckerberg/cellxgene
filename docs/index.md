_cellxgene_ is an interactive data explorer for single-cell transcriptomics data. Whether you need to visualize one thousand cells or one million, _cellxgene_ helps you gain insight into your single-cell data.

## features

#### flexible selections, coloring, and differential expression of your selected sets of cells
<img src="diffexp.gif" width="600"/>

#### single-gene analyses (e.g. expression analysis)
<img src="customGene.gif" width="600" />

## quick start

To install _cellxgene_ you need Python 3.6+. We recommend [installing _cellxgene_ into a conda or virtual environment.](/faq.html#how-do-i-create-a-python-environment-for-cellxgene)

Install the package.
``` bash
pip install cellxgene
```

Download an example [anndata](https://anndata.readthedocs.io/en/latest/) file

``` bash
curl -o pbmc3k.h5ad https://raw.githubusercontent.com/chanzuckerberg/cellxgene/master/example-dataset/pbmc3k.h5ad
```

Launch _cellxgene_
``` bash
cellxgene launch pbmc3k.h5ad --open
```

To explore more datasets already formatted for _cellxgene_, see [Data](data) or
visit [Getting Started](getting-started) to learn more about formatting your own
data for _cellxgene_.

## getting help

We'd love to hear from you!

For questions, suggestions, or accolades, [join the `#cellxgene-users` channel on the CZI Science Slack](https://join-cellxgene-users.herokuapp.com/) and say "hi!".

For any errors, [report bugs on Github](https://github.com/chanzuckerberg/cellxgene/issues).
