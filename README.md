<img src="./docs/cellxgene-logo.png" width="300">

_an interactive explorer for single-cell transcriptomics data_

[![DOI](https://zenodo.org/badge/105615409.svg)](https://zenodo.org/badge/latestdoi/105615409) [![PyPI](https://img.shields.io/pypi/v/cellxgene)](https://pypi.org/project/cellxgene/) [![PyPI - Downloads](https://img.shields.io/pypi/dm/cellxgene)](https://pypistats.org/packages/cellxgene) [![GitHub last commit](https://img.shields.io/github/last-commit/chanzuckerberg/cellxgene)](https://github.com/chanzuckerberg/cellxgene/pulse)
[![Push Tests](https://github.com/chanzuckerberg/cellxgene/workflows/Push%20Tests/badge.svg)](https://github.com/chanzuckerberg/cellxgene/actions?query=workflow%3A%22Push+Tests%22)
[![Compatibility Tests](https://github.com/chanzuckerberg/cellxgene/workflows/Compatibility%20Tests/badge.svg)](https://github.com/chanzuckerberg/cellxgene/actions?query=workflow%3A%22Compatibility+Tests%22)
![Code Coverage](https://codecov.io/gh/chanzuckerberg/cellxgene/branch/main/graph/badge.svg)

cellxgene (pronounced "cell-by-gene") is an interactive data explorer for single-cell transcriptomics datasets, such as those coming from the [Human Cell Atlas](https://humancellatlas.org). Leveraging modern web development techniques to enable fast visualizations of at least 1 million cells, we hope to enable biologists and computational researchers to explore their data.

Whether you need to visualize one thousand cells or one million, cellxgene helps you gain insight into your single-cell data.

<img src="https://github.com/chanzuckerberg/cellxgene/raw/main/docs/images/crossfilter.gif" width="350" height="200" hspace="30"><img src="https://github.com/chanzuckerberg/cellxgene/raw/main/docs/images/category-breakdown.gif" width="350" height="200" hspace="30">

# Getting started

### The comprehensive guide to cellxgene

[The cellxgene documentation is your one-stop-shop for information about cellxgene](https://chanzuckerberg.github.io/cellxgene/)! You may be particularly interested in:

- Seeing [what cellxgene can do](https://chanzuckerberg.github.io/cellxgene/posts/gallery)
- Learning more about cellxgene [installation](https://chanzuckerberg.github.io/cellxgene/posts/install) and [usage](https://chanzuckerberg.github.io/cellxgene/posts/launch)
- [Preparing your own data](https://chanzuckerberg.github.io/cellxgene/posts/prepare) for use in cellxgene
- Checking out [our roadmap](https://chanzuckerberg.github.io/cellxgene/posts/roadmap) for future development
- [Contributing](https://chanzuckerberg.github.io/cellxgene/posts/contribute) to cellxgene

### Quick start

To install cellxgene you need Python 3.6+. We recommend [installing cellxgene into a conda or virtual environment.](https://chanzuckerberg.github.io/cellxgene/posts/install)

Install the package.

```bash
pip install cellxgene
```

Launch cellxgene with an example [anndata](https://anndata.readthedocs.io/en/latest/) file

```bash
cellxgene launch https://cellxgene-example-data.czi.technology/pbmc3k.h5ad
```

To explore more datasets already formatted for cellxgene, check out the [Demo data](https://chanzuckerberg.github.io/cellxgene/posts/demo-data) or
see [Preparing your data](https://chanzuckerberg.github.io/cellxgene/posts/prepare) to learn more about formatting your own
data for cellxgene.

### Supported browsers

cellxgene currently supports the following browsers:

- Google Chrome 61+
- Edge 15+
- Firefox 60+
- Safari 10.1+

Please [file an issue](https://github.com/chanzuckerberg/cellxgene/issues/new/choose) if you would like us to add support for an unsupported browser.

### Finding help

We'd love to hear from you!
For questions, suggestions, or accolades, [join the `#cellxgene-users` channel on the CZI Science Slack](https://join-cellxgene-users.herokuapp.com/) and say "hi!".

For any errors, [report bugs on Github](https://github.com/chanzuckerberg/cellxgene/issues).

# Developing with cellxgene

### Contributing

We warmly welcome contributions from the community! Please see our [contributing guide](https://chanzuckerberg.github.io/cellxgene/posts/contribute) and don't hesitate to open an issue or send a pull request to improve cellxgene.

This project adheres to the Contributor Covenant [code of conduct](https://github.com/chanzuckerberg/.github/blob/master/CODE_OF_CONDUCT.md). By participating, you are expected to uphold this code. Please report unacceptable behavior to opensource@chanzuckerberg.com.

### Reuse

This project was started with the sole goal of empowering the scientific community to explore and understand their data. As such, we encourage other scientific tool builders in academia or industry to adopt the patterns, tools, and code from this project, and reach out to us with ideas or questions. All code is freely available for reuse under the [MIT license](https://opensource.org/licenses/MIT).

### Security

If you believe you have found a security issue, we would appreciate notification. Please send email to <security@chanzuckerberg.com>.

# Inspiration

We've been heavily inspired by several other related single-cell visualization projects, including the [UCSC Cell Browswer](http://cells.ucsc.edu/), [Cytoscape](http://www.cytoscape.org/), [Xena](https://xena.ucsc.edu/), [ASAP](https://asap.epfl.ch/), [Gene Pattern](http://genepattern-notebook.org/), and many others. We hope to explore collaborations where useful as this community works together on improving interactive visualization for single-cell data.

We were inspired by Mike Bostock and the [crossfilter](https://github.com/crossfilter) team for the design of our filtering implementation.

We have been working closely with the [scanpy](https://github.com/theislab/scanpy) team to integrate with their awesome analysis tools. Special thanks to Alex Wolf, Fabian Theis, and the rest of the team for their help during development and for providing an example dataset.

We are eager to explore integrations with other computational backends such as [Seurat](https://github.com/satijalab/seurat) or [Bioconductor](https://github.com/Bioconductor)
