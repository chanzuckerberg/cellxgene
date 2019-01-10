# cellxgene

> an interactive explorer for single-cell transcriptomics data

`cellxgene` is an interactive data explorer for single-cell transcriptomics datasets, such as those coming from the [Human Cell Atlas](https://humancellatlas.org). Leveraging modern web development techniques to enable fast visualizations of at least 1 million cells, we hope to enable biologists and computational researchers to explore their data, and to demonstrate general, scalable, and reusable patterns for scientific data visualization.

<img src="https://raw.githubusercontent.com/chanzuckerberg/cellxgene/master/docs/cellxgene-demo-1.gif" width="200" height="200" hspace="30"><img src="https://raw.githubusercontent.com/chanzuckerberg/cellxgene/master/docs/cellxgene-demo-2.gif" width="200" height="200" hspace="30"><img src="https://raw.githubusercontent.com/chanzuckerberg/cellxgene/master/docs/cellxgene-demo-3.gif" width="200" height="200" hspace="30">

## getting started

You'll need **python 3.6** and **Google Chrome**. (_Warning_: Python 3.7 is **not** supported at this time)
The web UI is tested on OSX and Windows using Chrome, and the python CLI is tested on OSX and Ubuntu (via WSL/Windows). It should work on other platforms, but if you run into trouble let us know (see [help](#help-and-contact) below).

To install run

```
pip install cellxgene
```

To start exploring a dataset call

```
cellxgene launch dataset.h5ad --open
```

If you want an example dataset download [this file](https://github.com/chanzuckerberg/cellxgene/raw/master/example-dataset/pbmc3k.h5ad) and then call

```
cellxgene launch pbmc3k.h5ad --open
```

You should see your web browser open with the following

<img width="450" src="https://raw.githubusercontent.com/chanzuckerberg/cellxgene/master/docs/cellxgene-opening-screenshot.png" pad="50px">

**Note**: automatic opening of the browser with the `--open` flag only works on OS X, on other platforms you'll need to directly point to the provided link in your browser.

There are several options available, such as:

- `--layout` to specify the layout as `tsne` or `umap`
- `--title` to show a title on the explorer
- `--open` to automatically open the web browser after launching (OS X only)

To see all options call

```
cellxgene launch --help
```

There is an additional subcommand called `cellxgene prepare` that takes an existing dataset in one of several formats and applies minimal preprocessing and reformatting so that `launch` can use it (see [the next section](##data-formatting) for more info on `prepare`).

## data formatting

### assumptions

The `launch` command assumes that the data is stored in the `.h5ad` format from the [`anndata`](https://anndata.readthedocs.io/en/latest/index.html) library. It also assumes that certain computations have already been performed. Briefly, the `.h5ad` format wraps a two-dimensional `ndarray` and stores additional metadata as "annotations" for either observations (referred to as `obs` and `obsm`) or variables (`var` and `varm`). `cellxgene launch` makes the following assumptions about your data (we recommend loading and inspecting your data using `scanpy` to validate these assumptions)

- an `obs` field has a unique identifier for every cell (you can specify which field to use with the `--obs-names` option, by default it will use the value of `data.obs_names`)
- a `var` field has a unique identifier for every gene (you can specify which field to use with the `--var-names` option, by default it will use the value of `data.var_names`)
- an `obsm` field contains the two-dimensional coordinates for the layout that you want to render (e.g. `X_tsne` for the `tsne` layout or `X_umap` for the `umap` layout)
- any additional `obs` fields will be rendered as per-cell continuous or categorical metadata by the app (e.g. `louvain` cluster assignments)

### prepare

The `prepare` command is included to help you format your data. It uses `scanpy` under the hood. This is especially useful if you are starting with raw unanalyzed data and are unfamiliar with `scanpy`.

To prepare from an existing `.h5ad` file use

```
cellxgene prepare dataset.h5ad --output=dataset-processed.h5ad
```

This will load the input data, perform PCA and nearest neighbor calculations, compute `umap` and `tsne` layouts and `louvain` cluster assignments, and save the results in a new file called `dataset-processed.h5ad` that can be loaded using `cellxgene launch`. Data can be loaded from several formats, including `.h5ad` `.loom` and a `10-Genomics-formatted` `mtx` directory. Several options are available, including running one of the preprocessing `recipes` included with `scanpy`, which include steps like cell filtering and gene selection. To learn more about the `recipes` please see the `scanpy` [documentation](https://github.com/theislab/scanpy/blob/master/scanpy/preprocessing/recipes.py).

Depending on the options chosen, `prepare` can take a long time to run (a few minutes for datasets with 10-100k cells, up to an hour or more for datasets with >100k cells). If you want `prepare` to run faster we recommend using the `sparse` option and only computing the layout for `umap`, using a call like this

```
cellxgene prepare dataset.h5ad --output=dataset-processed.h5ad --layout=umap --sparse
```

To see all options call

```
cellxgene prepare --help
```

**Note**: `cellxgene prepare` will only perform `louvain` clustering if you have the `python-igraph` and `louvain` packages installed. To make sure they are installed alongside `cellxgene` use

```
pip install cellxgene[louvain]
```

If the aforementioned optional package installation fails, you can also install these packages directly:

```
pip install python-igraph louvain>=0.6
```

## conda and virtual environments

If you use conda and want to create a conda environment for `cellxgene` you can use the following commands

```
conda create --yes -n cellxgene python=3.6
conda activate cellxgene
pip install cellxgene
```

Or you can create a virtual environment by using

```
ENV_NAME=cellxgene
python3.6 -m venv ${ENV_NAME}
source ${ENV_NAME}/bin/activate
pip install cellxgene
```

## docker

We have included a dockerfile to conveniently run cellxgene from docker.

1. Build the image `docker build . -t cellxgene`
2. Run the container and mount data `docker run -v "$PWD/example-dataset/:/data/" -p 5005:5005 cellxgene launch --host 0.0.0.0 data/pbmc3k.h5ad`
   - You will need to use --host 0.0.0.0 to have the container listen to incoming requests from the browser

## FAQ

<details>

<summary> questions about data formatting </summary>
  
<hr>
  
> Someone sent me a directory of `10X-Genomics` data with a `mtx` file and I've never used `scanpy`, can I use `cellxgene`?

Yep! This should only take a couple steps. We'll assume your data is in a folder called `data/` and you've successfully installed `cellxgene` with the `louvain` packages as described above. Just run

```
cellxgene prepare data/ --output=data-processed.h5ad --layout=umap
```

Depending on the size of the dataset, this may take some time. Once it's done, call

```
cellxgene launch data-processed.h5ad --layout=umap --open
```

And your web browser should open with an interactive view of your data.

<hr>

> In my `prepare` command I received the following error `Warning: louvain module is not installed, no clusters will be calculated. To fix this please install cellxgene with the optional feature louvain enabled`

Louvain clustering requires additional dependencies that are somewhat complex, so we don't include them by default. For now, you need to specify that you want these packages by using

```
pip install cellxgene[louvain]
```

<hr>

> I ran `prepare` and I'm getting results that look unexpected

You might want to try running one of the preprocessing recipes included with `scanpy` (read more about them [here](https://scanpy.readthedocs.io/en/latest/api/index.html#recipes)). You can specify this with the `--recipe` option, such as

```
cellxgene prepare data/ --output=data-processed.h5ad --recipe=zheng17
```

It should be easy to run `prepare` then call `cellxgene launch` a few times with different settings to explore different behaviors. We may explore adding other preprocessing options in the future.

<hr>

> I have extra metadata that I want to add to my dataset

Currently this is not supported directly, but you should be able to do this manually using `scanpy`. For example, this [notebook](https://github.com/falexwolf/fun-analyses/blob/master/tabula_muris/tabula_muris.ipynb) shows adding the contents of a `csv` file with metadata to an `anndata` object. For now, you could do this manually on your data in the same way and then save out the result before loading into `cellxgene`.

<hr>

> What part of the anndata objects does cellxgene pull in for visualization?

- `.obs` and `.var` annotations are use to extract metadata for filtering
- `.X` is used to display expression (histograms, scatterplot & colorscale) and to compute differential expression
- `.obsm` is used for layout

</details>

<details>
  
<summary> questions about installing and building </summary>

<hr>

> I tried to `pip install cellxgene` and got a weird error I don't understand

This may happen, especially as we work out bugs in our installation process! Please create a new [Github issue](https://github.com/chanzuckerberg/cellxgene/issues), explain what you did, and include all the error messages you saw. It'd also be super helpful if you call `pip freeze` and include the full output alongside your issue.

<hr>

> I'm following the developer instructions and get an error about "missing files and directories” when trying to build the client

This is likely because you do not have node and npm installed, we recommend using [nvm](https://github.com/creationix/nvm) if you're new to using these tools.

</details>

<details>

<summary> questions about algorithms </summary>

<hr>

> How are you computing and sorting differential expression results?

Currently we use a [Welch's _t_-test](https://en.wikipedia.org/wiki/Welch%27s_t-test) implementation including the same variance overestimation correction as used in `scanpy`. We sort the `tscore` to identify the top N genes, and then filter to remove any that fall below a cutoff log fold change value, which can help remove spurious test results. The default threshold is `0.01` and can be changed using the option `--diffexp-lfc-cutoff`. We can explore adding support for other test types in the future.

</details>

## developer guide

This project has made a few key design choices

- The front-end is built with [`regl`](https://github.com/regl-project/regl) (a webgl library), [`react`](https://reactjs.org/), [`redux`](https://redux.js.org/), [`d3`](https://github.com/d3/d3), and [`blueprint`](https://blueprintjs.com/docs/#core) to handle rendering large numbers of cells with lots of complex interactivity
- The app is designed with a client-server model that can support a range of existing analysis packages for backend computational tasks (currently built for [scanpy](https://github.com/theislab/scanpy))
- The client uses fast cross-filtering to handle selections and comparisons across subsets of data

Depending on your background and interests, you might want to contribute to the frontend, or backend, or both!

If you are interested in working on `cellxgene` development, we recommend cloning the project from Gitub. First you'll need the following installed on your machine

- python 3.6
- node and npm (we recommend using [nvm](https://github.com/creationix/nvm) if this is your first time with node)

Then clone the project

```
git clone https://github.com/chanzuckerberg/cellxgene.git
```

Build the client web assets by calling this from inside the `cellxgene` folder

```
./bin/build-client
```

Install all requirements (we recommend doing this inside a virtual environment)

```
pip install -e .
```

You can start the app while developing either by calling `cellxgene` or by calling `python -m server`. We recommend using the `--debug` flag to see more output, which you can include when reporting bugs.

If you have any questions about developing or contributing, come hang out with us by joining the [CZI Science Slack](https://join-cellxgene-users.herokuapp.com/) and posting in the `#cellxgene-dev` channel.

## development roadmap

`cellxgene` is still very much in development, and we've love to include the community as we plan new features to work on. We are thinking about working on the following features over the next 3-12 months. If you are interested in updates, want to give feedback, want to contribute, or have ideas about other features we should work on, please [contact us](#help-and-contact)

- **Visualizaling spatial metadata** Image-based transcriptomics methods also generate large cell by gene matrices, alongside rich metadata about spatial location; we would like to render this information in `cellxgene`
- **Visualizing trajectories** Trajectory analyses infer progression along some ordering or pseudotime; we would like `cellxgene` to render the results of these analyses when they have been performed
- **Deploy to web** Many projects release public data browser websites alongside their publicatons; we would like to make it easy for anyone to deploy `cellxgene` to a custom URL with their own dataset that they own and operate
- **HCA Integration** The [Human Cell Atlas](https://humancellatlas.org) is generating a large corpus of single-cell expression data and will make it available through the Data Coordination Platform; we would like `cellxgene` to be one of several different portals for browsing these data

## contributing

We warmly welcome contributions from the community! Please submit any bug reports and feature requests through [Github issues](https://github.com/chanzuckerberg/cellxgene/issues). Please submit any direct contributions by forking the repository, creating a branch, and submitting a Pull Request. It'd be great for PRs to include test cases and documentation updates where relevant, though we know the core test suite is itself still a work in progress. And all code contributions and dependencies must be compatible with the project's open-source license (MIT). If you have any questions about this stuff, just ask!

## inspiration and collaboration

We've been heavily inspired by several other related single-cell visualization projects, including the [UCSC Cell Browswer](http://cells.ucsc.edu/), [Cytoscape](http://www.cytoscape.org/), [Xena](https://xena.ucsc.edu/), [ASAP](https://asap.epfl.ch/), [Gene Pattern](http://genepattern-notebook.org/), and many others. We hope to explore collaborations where useful as this community works together on improving interactive visualization for single-cell data.

We were inspired by Mike Bostock and the [crossfilter](https://github.com/crossfilter) team for the design of our filtering implementation.

We have been working closely with the [`scanpy`](https://github.com/theislab/scanpy) team to integrate with their awesome analysis tools. Special thanks to Alex Wolf, Fabian Theis, and the rest of the team for their help during development and for providing an example dataset.

We are eager to explore integrations with other computational backends such as [`Seurat`](https://github.com/satijalab/seurat) or [`Bioconductor`](https://github.com/Bioconductor)

## help and contact

Have questions, suggestions, or comments? You can come hang out with us by joining the [CZI Science Slack](https://join-cellxgene-users.herokuapp.com/) and posting in the `#cellxgene-users` channel. As mentioned above, please submit any feature requests or bugs as [Github issues](https://github.com/chanzuckerberg/cellxgene/issues). We'd love to hear from you!

## reuse

This project was started with the sole goal of empowering the scientific community to explore and understand their data. As such, we encourage other scientific tool builders in academia or industry to adopt the patterns, tools, and code from this project, and reach out to us with ideas or questions. All code is freely available for reuse under the [MIT license](https://opensource.org/licenses/MIT).
