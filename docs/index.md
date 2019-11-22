---
title: Index
subtitle: Index
layout: default
---
# Quick start

Whether you need to visualize one thousand cells or one million, cellxgene helps you gain insight into your single-cell data.

To install cellxgene you need Python 3.6+. We recommend [installing cellxgene into a conda or virtual environment.](/faq.html#how-do-i-create-a-python-environment-for-cellxgene)

Install the package.
``` bash
pip install cellxgene
```

Download an example [anndata](https://anndata.readthedocs.io/en/latest/) file

``` bash
curl -o tabula-muris.h5ad https://cellxgene-example-data.czi.technology/tabula-muris.h5ad.zip
unzip tabula-muris.h5ad.zip
```

Launch cellxgene
``` bash
cellxgene launch tabula-muris.h5ad --open
```

To explore more datasets already formatted for cellxgene, check out the [Demo data](demo-data) or
see [Preparing your data](prepare) to learn more about formatting your own
data for cellxgene.

# Getting help

We'd love to hear from you!

For questions, suggestions, or accolades, [join the `#cellxgene-users` channel on the CZI Science Slack](https://join-cellxgene-users.herokuapp.com/) and say "hi!".

For any errors, [report bugs on Github](https://github.com/chanzuckerberg/cellxgene/issues).
