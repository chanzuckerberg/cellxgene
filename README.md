<img src="./docs/cellxgene-logo.png" width="300">

_an interactive explorer for single-cell transcriptomics data_

[![DOI](https://zenodo.org/badge/105615409.svg)](https://zenodo.org/badge/latestdoi/105615409) [![PyPI](https://img.shields.io/pypi/v/cellxgene)](https://pypi.org/project/cellxgene/) [![PyPI - Downloads](https://img.shields.io/pypi/dm/cellxgene)](https://pypistats.org/packages/cellxgene) [![GitHub last commit](https://img.shields.io/github/last-commit/chanzuckerberg/cellxgene)](https://github.com/chanzuckerberg/cellxgene/pulse)
[![Push Tests](https://github.com/chanzuckerberg/cellxgene/workflows/Push%20Tests/badge.svg)](https://github.com/chanzuckerberg/cellxgene/actions?query=workflow%3A%22Push+Tests%22)
[![Compatibility Tests](https://github.com/chanzuckerberg/cellxgene/workflows/Compatibility%20Tests/badge.svg)](https://github.com/chanzuckerberg/cellxgene/actions?query=workflow%3A%22Compatibility+Tests%22)
![Code Coverage](https://codecov.io/gh/chanzuckerberg/cellxgene/branch/main/graph/badge.svg)

cellxgene Desktop (pronounced "cell-by-gene") is an interactive data explorer for single-cell datasets, such as those coming from the [Human Cell Atlas](https://humancellatlas.org). Leveraging modern web development techniques to enable fast visualizations of at least 1 million cells, we hope to enable biologists and computational researchers to explore their data.

Whether you need to visualize one thousand cells or one million, cellxgene Desktop helps you gain insight into your single-cell data.

<img src="https://github.com/chanzuckerberg/cellxgene/raw/main/docs/images/crossfilter.gif" width="350" height="200" hspace="30"><img src="https://github.com/chanzuckerberg/cellxgene/raw/main/docs/images/category-breakdown.gif" width="350" height="200" hspace="30">

# Getting started

### The comprehensive guide to cellxgene Desktop

[The cellxgene documentation is your one-stop-shop for information about cellxgene Desktop](https://github.com/chanzuckerberg/cellxgene-documentation/blob/main/README.md)! You may be particularly interested in:

- Seeing [what cellxgene Desktop can do](https://github.com/chanzuckerberg/cellxgene-documentation/blob/main/explore-data/explorer-tutorials.md)
- Learning more about cellxgene [installation](https://github.com/chanzuckerberg/cellxgene-documentation/blob/main/desktop/install.md) and [usage](https://github.com/chanzuckerberg/cellxgene-documentation/blob/main/desktop/quick-start.md#quick-start-1)
- [Preparing your own data](https://github.com/chanzuckerberg/cellxgene-documentation/blob/main/desktop/data-reqs.md) for use in cellxgene Desktop
- Checking out [our roadmap](https://github.com/chanzuckerberg/cellxgene-documentation/blob/main/roadmap.md) for future development
- [Contributing](https://github.com/chanzuckerberg/cellxgene-documentation/blob/main/contribute.md) to cellxgene Desktop

### Quick start

To install cellxgene Desktop you need Python 3.6+. We recommend [installing cellxgene Desktop into a conda or virtual environment.](https://github.com/chanzuckerberg/cellxgene-documentation/blob/main/desktop/install.md)

Install the package.

```bash
pip install cellxgene
```

Launch cellxgene Desktop with an example [anndata](https://anndata.readthedocs.io/en/latest/) file

```bash
cellxgene launch https://cellxgene-example-data.czi.technology/pbmc3k.h5ad
```

To explore more datasets already formatted for cellxgene Desktop, check out the [Demo data](https://github.com/chanzuckerberg/cellxgene-documentation/blob/main/desktop/quick-start.md#example-datasets) or
see [Preparing your data](https://github.com/chanzuckerberg/cellxgene-documentation/blob/main/desktop/data-reqs.md) to learn more about formatting your own
data for cellxgene Desktop.

### Supported browsers

cellxgene Desktop currently supports the following browsers:

- Google Chrome 61+
- Edge 15+
- Firefox 60+

Please [file an issue](https://github.com/chanzuckerberg/cellxgene/issues/new/choose) if you would like us to add support for an unsupported browser.

### Finding help

We'd love to hear from you!
For questions, suggestions, or accolades, [join the `#cellxgene-users` channel on the CZI Science Slack](https://join-cellxgene-users.herokuapp.com/) and say "hi!".

For any errors, [report bugs on Github](https://github.com/chanzuckerberg/cellxgene/issues).

# Developing with cellxgene Desktop

### Contributing

We warmly welcome contributions from the community! Please see our [contributing guide](https://github.com/chanzuckerberg/cellxgene-documentation/blob/main/contribute.md) and don't hesitate to open an issue or send a pull request to improve cellxgene Desktop. Please see the [dev_docs](https://github.com/chanzuckerberg/cellxgene/tree/main/dev_docs) for pull request suggestions, unit test details, local documentation preview, and other development specifics. 

This project adheres to the Contributor Covenant [code of conduct](https://github.com/chanzuckerberg/.github/blob/master/CODE_OF_CONDUCT.md). By participating, you are expected to uphold this code. Please report unacceptable behavior to opensource@chanzuckerberg.com.

### Reuse

This project was started with the sole goal of empowering the scientific community to explore and understand their data. 
As such, we encourage other scientific tool builders in academia or industry to adopt the patterns, tools, and code from 
this project. All code is freely available for reuse under the [MIT license](https://opensource.org/licenses/MIT).


Before extending cellxgene, we encourage you to reach out to us with ideas or questions. It might be possible that an 
extension could be directly contributed, which would make it available for a wider audience, or that it's on our 
[roadmap](https://github.com/chanzuckerberg/cellxgene-documentation/blob/main/roadmap.md) and under active development. 

See the [cellxgene extensions](https://github.com/chanzuckerberg/cellxgene-documentation/blob/main/community-extensions.md) section of our documentation for examples of community use and cellxgene extensions. 

### Security

If you believe you have found a security issue, we would appreciate notification. Please send email to <security@chanzuckerberg.com>.

# Inspiration

We've been heavily inspired by several other related single-cell visualization projects, including the [UCSC Cell Browser](http://cells.ucsc.edu/), [Cytoscape](http://www.cytoscape.org/), [Xena](https://xena.ucsc.edu/), [ASAP](https://asap.epfl.ch/), [GenePattern](http://genepattern-notebook.org/), and many others. We hope to explore collaborations where useful as this community works together on improving interactive visualization for single-cell data.

We were inspired by Mike Bostock and the [crossfilter](https://github.com/crossfilter) team for the design of our filtering implementation.

We have been working closely with the [scanpy](https://github.com/theislab/scanpy) team to integrate with their awesome analysis tools. Special thanks to Alex Wolf, Fabian Theis, and the rest of the team for their help during development and for providing an example dataset.

We are eager to explore integrations with other computational backends such as [Seurat](https://github.com/satijalab/seurat) or [Bioconductor](https://github.com/Bioconductor)
