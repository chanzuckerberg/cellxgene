_cellxgene_ is an interactive data explorer for single-cell transcriptomics data designed to handle large datasets (1 million cells or more) and integrate with your favorite analysis tools

## features

- Flexible selection and coloring of cells
- Differential expression of arbitrary sets of cells
- Single-gene analyses

## getting started

Install the package
``` bash
pip install cellxgene
```

Download an example file

``` bash
curl -o pbmc3k.h5ad https://raw.githubusercontent.com/chanzuckerberg/cellxgene/master/example-dataset/pbmc3k.h5ad
```

Launch _cellxgene_
``` bash
cellxgene launch pbmc3k.h5ad
```

## getting help

We'd love to hear from you!

For questions, suggestions, or accolades, [join the `#cellxgene-users` channel on the CZI Science Slack](https://join-cziscience-slack.herokuapp.com/) and say "hi!".

For any errors, [report bugs on Github](https://github.com/chanzuckerberg/cellxgene/issues).
