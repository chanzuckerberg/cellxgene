<img src="./docs/cellxgene-logo.svg" width="300">

_an interactive explorer for single-cell transcriptomics data_

[![DOI](https://zenodo.org/badge/105615409.svg)](https://zenodo.org/badge/latestdoi/105615409) [![PyPI](https://img.shields.io/pypi/v/cellxgene)](https://pypi.org/project/cellxgene/) [![PyPI - Downloads](https://img.shields.io/pypi/dm/cellxgene)](https://pypistats.org/packages/cellxgene) [![GitHub last commit](https://img.shields.io/github/last-commit/chanzuckerberg/cellxgene)](https://github.com/chanzuckerberg/cellxgene/pulse)

cellxgene (pronounced "cell-by-gene") is an interactive data explorer for single-cell transcriptomics datasets, such as those coming from the [Human Cell Atlas](https://humancellatlas.org). Leveraging modern web development techniques to enable fast visualizations of at least 1 million cells, we hope to enable biologists and computational researchers to explore their data.

<img src="https://raw.githubusercontent.com/chanzuckerberg/cellxgene/master/docs/cellxgene-demo-1.gif" width="200" height="200" hspace="30"><img src="https://raw.githubusercontent.com/chanzuckerberg/cellxgene/master/docs/cellxgene-demo-2.gif" width="200" height="200" hspace="30"><img src="https://raw.githubusercontent.com/chanzuckerberg/cellxgene/master/docs/cellxgene-demo-3.gif" width="200" height="200" hspace="30">

- Want to install and use cellxgene? Visit the [cellxgene docs](https://chanzuckerberg.github.io/cellxgene/).
- Want to see where we are going? Check out [our roadmap](https://github.com/chanzuckerberg/cellxgene/blob/master/ROADMAP.md).
- Want to contribute? See our [contributors guide](https://github.com/chanzuckerberg/cellxgene/blob/master/CONTRIBUTING.md).

## quick start

To install cellxgene you need Python 3.6+. We recommend [installing cellxgene into a conda or virtual environment.](https://chanzuckerberg.github.io/cellxgene/faq.html#how-do-i-create-a-python-36-environment-for-cellxgene)

Install the package.

```bash
pip install cellxgene
```

Download an example [anndata](https://anndata.readthedocs.io/en/latest/) file

```bash
curl -O https://cellxgene-example-data.czi.technology/pbmc3k.h5ad.zip
unzip pbmc3k.h5ad
```

Launch cellxgene

```bash
cellxgene launch pbmc3k.h5ad --open
```

To learn more about what you can do with cellxgene, see the [Getting Started](https://chanzuckerberg.github.io/cellxgene/getting-started.html) guide.

## get in touch

Have questions, suggestions, or comments? You can come hang out with us by joining the [CZI Science Slack](https://join-cellxgene-users.herokuapp.com/) and posting in the `#cellxgene-users` channel. Have feature requests or bugs? Please submit these as [Github issues](https://github.com/chanzuckerberg/cellxgene/issues). We'd love to hear from you!

## contributing

We warmly welcome contributions from the community! Please see our [contributing guide](https://github.com/chanzuckerberg/cellxgene/blob/master/CONTRIBUTING.md) and don't hesitate to open an issue or send a pull request to improve cellxgene.

This project adheres to the Contributor Covenant [code of conduct](https://github.com/chanzuckerberg/.github/blob/master/CODE_OF_CONDUCT.md). By participating, you are expected to uphold this code. Please report unacceptable behavior to opensource@chanzuckerberg.com.

## core team

The current core team:

- Colin Megill, frontend & product design
- Bruce Martin, software engineer
- Sidney Bell, computational biologist
- Lia Prins, designer
- Severiano Badajoz, software engineer

We would also like to gratefully acknowledge contributions from past core team members:

- Charlotte Weaver, software engineer

## where we are going

Our goal is to enable teams of computational and experimental
biologists to collaboratively gain insight into their single-cell RNA-seq data.

There are 4 key features we plan to implement in the near term.

- Click install and launch
- Manual annotation workflows
- Toggle embeddings
- Gene information

For more detail on these features and where we are going, see [our roadmap](https://github.com/chanzuckerberg/cellxgene/blob/master/ROADMAP.md).

## inspiration

We've been heavily inspired by several other related single-cell visualization projects, including the [UCSC Cell Browswer](http://cells.ucsc.edu/), [Cytoscape](http://www.cytoscape.org/), [Xena](https://xena.ucsc.edu/), [ASAP](https://asap.epfl.ch/), [Gene Pattern](http://genepattern-notebook.org/), and many others. We hope to explore collaborations where useful as this community works together on improving interactive visualization for single-cell data.

We were inspired by Mike Bostock and the [crossfilter](https://github.com/crossfilter) team for the design of our filtering implementation.

We have been working closely with the [scanpy](https://github.com/theislab/scanpy) team to integrate with their awesome analysis tools. Special thanks to Alex Wolf, Fabian Theis, and the rest of the team for their help during development and for providing an example dataset.

We are eager to explore integrations with other computational backends such as [Seurat](https://github.com/satijalab/seurat) or [Bioconductor](https://github.com/Bioconductor)

## reuse

This project was started with the sole goal of empowering the scientific community to explore and understand their data. As such, we encourage other scientific tool builders in academia or industry to adopt the patterns, tools, and code from this project, and reach out to us with ideas or questions. All code is freely available for reuse under the [MIT license](https://opensource.org/licenses/MIT).

## security

If you believe you have found a security issue, we would appreciate notification. Please send email to <security@chanzuckerberg.com>.
