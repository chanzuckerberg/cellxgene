# cellxgene

> an interactive performant explorer for single cell transcriptomics data

`cellxgene` is an interactive data explorer that leverages modern web development techniques to visualize large single-cell transcriptomics datasets, such as the data coming from large consortia projects like the Human Cell Atlas. We want it to enable biologists and computational researchers to explore their data, and to demonstrate general, scalable, reusable patterns for visualizing large scientific datasets.

<img width="500" src="./docs/cellxgene-demo.gif" pad="50px">

## features

- _Visualization at scale_  Built with [WebGL](https://www.khronos.org/webgl/), [React](https://reactjs.org/), & [Redux](https://redux.js.org/) to enable performant exploration of at least 1 million cells.
- _Interactive exploration_  Select, cross-filter, and compare subsets of data with performant indexing and data handling.
- _Flexible API_  The cellxgene client-server model is designed to support a range of existing analysis packages for backend computational tasks (currently built for [scanpy](https://github.com/theislab/scanpy)), integrated with client-side visualization via a [REST API](https://restfulapi.net/).

## getting started

You'll need **python 3.6** or later and **Google Chrome**. Currently tested on OSX or Windows (via WSL using Ubuntu). It should work on other platforms, if you run into trouble let us know (see [help](#help-and-contact) below).

To install call

```
pip install cellxgene
```

To start the tool and see the help information call

```
cellxgene --help
```

There are two primary subcommands: `launch` and `prepare`:

- The `launch` command loads a dataset and starts the interactive explorer in your web browser
- The optional `prepare` command takes an existing dataset in one of several formats and applies minimal preprocessing and reformatting so that `launch` can use it

To explore a dataset call

```
cellxgene launch dataset.h5ad --open
```

If you want an example dataset download [this file](https://github.com/chanzuckerberg/cellxgene/raw/master/example-dataset/pbmc3k.h5ad) and then call

```
cellxgene launch pbmc3k.h5ad --open
```
You should see your web browser open with the following

<img width="350" height="218" src="./docs/cellxgene-screenshot.png" pad="50px">

There are several options available, such as:

- `--layout` to specify the layout as  `tsne` or `umap`
- `--title` to show a title on the explorer
- `--open` to automatically open the web browser after launching

To see all options call

```
cellxgene launch --help
```

## data formatting

### assumptions

The `launch` command assumes that the data is in the `.h5ad` `anndata` format (read more about [anndata](https://anndata.readthedocs.io/en/latest/index.html)). It makes the following assumptions about the structure of the data:

- an `obsm` field contains the coordinates for the layout that you want to render (e.g. `X_tsne` for the `tsne` layout or `X_umap` for the `umap` layout)
- an `obs`  annotation has a unique identifier for every cell (you can specify this using the `--obs-names` option, by default it will use the index)
- an `var` annotation has a unique identifier for every gene (you can specify this using the `--var-names` option, by default it will use the index)
- any additional  `obs` annotations will be rendered as per-cell metadata by the app (e.g. `louvain` cluster assignments) 

### prepare

The `prepare` command is included to help you format your data using `scanpy`. This is especially useful if you are unfamiliar with `scanpy` or have a data in a different format. 

To prepare from an existing `.h5ad` file use

```
cellxgene prepare dataset.h5ad --output=dataset-processed.h5ad
```

This will load the input data, perform PCA and nearest neighbors, compute `umap` and `tsne` layouts and `louvain` cluster assignments, and save the results in a new file called `dataset-processed.h5ad` that can be loaded using `cellxgene launch` . Several options are available, including running the preprocessing `recipes` included with `scanpy`. To see all options call

```
cellxgene prepare --help
```

_Note_: `cellxgene prepare` will only perform `louvain` clustering if you have the `python-igraph` and `louvain` packages installed. To make sure they are installed with `cellxgene` install using

```
pip install cellxgene[louvain]
```

## developer guide

If you are interested in working on `cellxgene` development, we recommend cloning the project from Gitub. First you'll need the following installed on your machine

- python 3.6
- node and npm (we recommend using [nvm](https://github.com/creationix/nvm) if this is your first time with node)

Then clone the project

```
git clone https://github.com/chanzuckerberg/cellxgene.git
```

Build the client assets by calling this from inside the `cellxgene` folder

```
./bin/build-client
```

Install all requirements (we recommend using virtual environments)

```
pip install -e .
```

You can start the app while developing either by calling `cellxgene` or by calling `python -m server`.

## faq

> I'm following the developer instructions and get an error about "missing files and directories” when trying to build the client

This is likely because you do not have node and npm installed, we recommend using [nvm](https://github.com/creationix/nvm)

## development roadmap

`cellxgene` is still very much in development, and we've love to include the community as we plan new features to work on. We are thinking about working on the following features over the next 3-12 months. If you are interested in updates, want to give feedback, want to contribute, or have ideas about other features we should work on, please [contact us](# Help/Contact) 

- _Visualizaling spatial metadata_  Image-based transcriptomics methods also generate large cell by gene matrices, alongside rich metadata about spatial location; we would like to render this information in `cellxgene`.
- _Visualizing trajectories_  Trajectory analyses infer progression along some ordering or pseudotime; we would like `cellxgene ` to render the results of these analyses when they have been performed
- _Deploy to web_  Many projects release public data browser websites  alongside their publicatons; we would like to make it easy for anyone to deploy `cellxgene` to a custom URL with their own dataset that they own and operate.
- _HCA Integration_  The [Human Cell Atlas](https://humancellatlas.org) is generating a large corupus of single-cell expression data and will make it available through the Data Coordination Platform; we would like `cellxgene` to be one of several different portals for browsing these data

## contributing

We warmly welcome contributions from the community! Please submit any bug reports and feature requests through [Github issues](https://github.com/chanzuckerberg/cellxgene/issues). Please submit any direct contributions by forking the repository, creating a branch (optional), and submitting a Pull Request.

## inspiration and collaboration

We’ve been heavily inspired by several other related single-cell visualization projects, including the [UCSC Cell Browswer](http://cells.ucsc.edu/), [Cytoscape](http://www.cytoscape.org/), [Xena](https://xena.ucsc.edu/), [ASAP](https://asap.epfl.ch/), [Gene Pattern](http://genepattern-notebook.org/), and many others. We hope to explore collaborations where useful as this community works together on improving interactive visualization for single-cell data.

We were inspired by Mike Bostock and the [crossfilter](https://github.com/crossfilter) team for the design of our filtering implementation.

We have been working closely with the [`scanpy`](https://github.com/theislab/scanpy) team to integrate with their awesome analysis tools. Special thanks to Alex Wolf, Fabian Theis, and the rest of the team for their help during development and for providing an example dataset. 

We are eager to explore integrations with other computational backends such as [`Seurat`](https://github.com/satijalab/seurat) or [`Bioconductor`](https://github.com/Bioconductor)

## help and contact

Have questions, suggestions, or comments? You can come hang out with us by joining the [CZI Science Slack](https://cziscience.slack.com/messages/CCTA8DF1T) and posting in the `#cellxgene` channel. As mentioned above, please submit any feature requests or bugs as [Github issues](https://github.com/chanzuckerberg/cellxgene/issues). We'd love to hear from you!

## reuse

This project was started with the sole goal of empowering the scientific community to explore and understand their data. As such, we encourage other scientific tool builders to adopt the patterns, tools, and code from this project, and reach out to us with ideas or questions. All code is freely available for reuse under the [MIT license](https://opensource.org/licenses/MIT).
