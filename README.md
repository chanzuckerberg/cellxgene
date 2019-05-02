# cellxgene

> an interactive explorer for single-cell transcriptomics data

_cellxgene_ (pronounced "sell-by-jean") is an interactive data explorer for single-cell transcriptomics datasets, such as those coming from the [Human Cell Atlas](https://humancellatlas.org). Leveraging modern web development techniques to enable fast visualizations of at least 1 million cells, we hope to enable biologists and computational researchers to explore their data, and to demonstrate general, scalable, and reusable patterns for scientific data visualization.

<img src="https://raw.githubusercontent.com/chanzuckerberg/cellxgene/master/docs/cellxgene-demo-1.gif" width="200" height="200" hspace="30"><img src="https://raw.githubusercontent.com/chanzuckerberg/cellxgene/master/docs/cellxgene-demo-2.gif" width="200" height="200" hspace="30"><img src="https://raw.githubusercontent.com/chanzuckerberg/cellxgene/master/docs/cellxgene-demo-3.gif" width="200" height="200" hspace="30">

- Want to install and use cellxgene? Visit the [cellxgene docs](https://chanzuckerberg.github.io/cellxgene/).
- Want to see where we are going? Check out [our roadmap](ROADMAP.md).
- Want to contribute? See our [contributors guide](#Contributing)

## quick start

To install _cellxgene_ you need Python 3.6+. We recommend [installing _cellxgene_ into a conda or virtual environment.](https://chanzuckerberg.github.io/cellxgene/faq.html#how-do-i-create-a-python-36-environment-for-cellxgene)

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

To learn more about what you can do with _cellxgene_, see the [Getting Started](https://chanzuckerberg.github.io/cellxgene/getting-stared/) guide.

## get in touch

Have questions, suggestions, or comments? You can come hang out with us by joining the [CZI Science Slack](https://join-cellxgene-users.herokuapp.com/) and posting in the `#cellxgene-users` channel. As mentioned above, please submit any feature requests or bugs as [Github issues](https://github.com/chanzuckerberg/cellxgene/issues). We'd love to hear from you!

## where we are going

Our goal is to enable teams of computational and experimental
biologists to collaboratively gain insight into their single-cell RNA-seq data.

There are 4 key features we plan to implement in the near term.

- Click install and launch
- Manual annotation workflows
- Toggle embeddings
- Gene information

For more detail on these features and where we are going, see [our roadmap](ROADMAP.md).

## contributing

We warmly welcome contributions from the community! Please submit any bug reports and feature requests through [Github issues](https://github.com/chanzuckerberg/cellxgene/issues). Please submit any direct contributions by forking the repository, creating a branch, and submitting a Pull Request. It'd be great for PRs to include test cases and documentation updates where relevant, though we know the core test suite is itself still a work in progress. And all code contributions and dependencies must be compatible with the project's open-source license (MIT). If you have any questions about this stuff, just ask!

### developer guide

This project has made a few key design choices

- The front-end is built with [`regl`](https://github.com/regl-project/regl) (a webgl library), [`react`](https://reactjs.org/), [`redux`](https://redux.js.org/), [`d3`](https://github.com/d3/d3), and [`blueprint`](https://blueprintjs.com/docs/#core) to handle rendering large numbers of cells with lots of complex interactivity
- The app is designed with a client-server model that can support a range of existing analysis packages for backend computational tasks (currently built for [scanpy](https://github.com/theislab/scanpy))
- The client uses fast cross-filtering to handle selections and comparisons across subsets of data

Depending on your background and interests, you might want to contribute to the frontend, or backend, or both!

If you are interested in working on `cellxgene` development, we recommend cloning the project from Gitub. First you'll need the following installed on your machine

- python 3.6+
- node and npm (we recommend using [nvm](https://github.com/creationix/nvm) if this is your first time with node)

Then clone the project

```
git clone https://github.com/chanzuckerberg/cellxgene.git
```

Build the client web assets by calling `make` from inside the `cellxgene` folder

```
make
```

Install all requirements (we recommend doing this inside a virtual environment)

```
pip install -e .
```

You can start the app while developing either by calling `cellxgene` or by calling `python -m server`. We recommend using the `--debug` flag to see more output, which you can include when reporting bugs.

If you have any questions about developing or contributing, come hang out with us by joining the [CZI Science Slack](https://join-cellxgene-users.herokuapp.com/) and posting in the `#cellxgene-dev` channel.

## inspiration

We've been heavily inspired by several other related single-cell visualization projects, including the [UCSC Cell Browswer](http://cells.ucsc.edu/), [Cytoscape](http://www.cytoscape.org/), [Xena](https://xena.ucsc.edu/), [ASAP](https://asap.epfl.ch/), [Gene Pattern](http://genepattern-notebook.org/), and many others. We hope to explore collaborations where useful as this community works together on improving interactive visualization for single-cell data.

We were inspired by Mike Bostock and the [crossfilter](https://github.com/crossfilter) team for the design of our filtering implementation.

We have been working closely with the [`scanpy`](https://github.com/theislab/scanpy) team to integrate with their awesome analysis tools. Special thanks to Alex Wolf, Fabian Theis, and the rest of the team for their help during development and for providing an example dataset.

We are eager to explore integrations with other computational backends such as [`Seurat`](https://github.com/satijalab/seurat) or [`Bioconductor`](https://github.com/Bioconductor)

## core team

- Colin Megill, frontend & product design
- Charlotte Weaver, software engineer
- Bruce Martin, software engineer
- Sidney Bell, computational biologist
- Justin Kiggins, product manager

## reuse

This project was started with the sole goal of empowering the scientific community to explore and understand their data. As such, we encourage other scientific tool builders in academia or industry to adopt the patterns, tools, and code from this project, and reach out to us with ideas or questions. All code is freely available for reuse under the [MIT license](https://opensource.org/licenses/MIT).
