## getting started

You'll need **python 3.6+** and **Google Chrome**.
The web UI is tested on OSX and Windows using Chrome, and the python CLI is tested on OSX and Ubuntu (via WSL/Windows). It should work on other platforms, but if you run into trouble let us know.

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

On Mac OS and Ubuntu, you should see your web browser open with the following

<img width="450" src="https://raw.githubusercontent.com/chanzuckerberg/cellxgene/master/docs/cellxgene-opening-screenshot.png" pad="50px">

**Note**: automatic opening of the browser with the `--open` flag only works on OS X, on other platforms you'll need to directly point to the provided link in your browser.

There are several options available, such as:

- `--embedding` to restrict available emdeddings in the UI, eg, `tsne`, `umap`, `diffmap`, `phate`, `draw_graph_fa`, or `draw_graph_fr`
- `--title` to show a title on the explorer
- `--open` to automatically open the web browser after launching (OS X only)

To see all options call

```
cellxgene launch --help
```

There is an additional subcommand called `cellxgene prepare` that takes an existing dataset in one of several formats and applies minimal preprocessing and reformatting so that `launch` can use it (see [the next section](##data-formatting) for more info on `prepare`).

## data formatting

### requirements

The `launch` command assumes that the data is stored in the `.h5ad` format from the [`anndata`](https://anndata.readthedocs.io/en/latest/index.html) library. It also assumes that certain computations have already been performed. Briefly, the `.h5ad` format wraps a two-dimensional `ndarray` and stores additional metadata as "annotations" for either observations (referred to as `obs` and `obsm`) or variables (`var` and `varm`). `cellxgene launch` makes the following assumptions about your data (we recommend loading and inspecting your data using `scanpy` to validate these assumptions)

- an `obs` field has a unique identifier for every cell (you can specify which field to use with the `--obs-names` option, by default it will use the value of `data.obs_names`)
- a `var` field has a unique identifier for every gene (you can specify which field to use with the `--var-names` option, by default it will use the value of `data.var_names`)
- an `obsm` field contains the two-dimensional coordinates for the embedding that you want to render (e.g. `X_umap` for the `umap` embedding)
- any additional `obs` fields will be rendered as per-cell continuous or categorical metadata by the app (e.g. `louvain` cluster assignments)

### prepare

The `prepare` command is included to help you format your data. It uses `scanpy` under the hood. This is especially useful if you are starting with raw unanalyzed data and are unfamiliar with `scanpy`.

To install `cellxgene prepare` alongside `cellxgene launch`

```
pip install cellxgene[prepare]
```

If the aforementioned optional package installation fails, you can also install these packages directly:

```
pip install scanpy>=1.3.7 python-igraph louvain>=0.6
```

To prepare from an existing `.h5ad` file use

```
cellxgene prepare dataset.h5ad --output=dataset-processed.h5ad
```

This will load the input data, perform PCA and nearest neighbor calculations, compute `umap` and `tsne` embeddings and `louvain` cluster assignments, and save the results in a new file called `dataset-processed.h5ad` that can be loaded using `cellxgene launch`. Data can be loaded from several formats, including `.h5ad` `.loom` and a `10-Genomics-formatted` `mtx` directory. Several options are available, including running one of the preprocessing `recipes` included with `scanpy`, which include steps like cell filtering and gene selection. To learn more about the `recipes` please see the `scanpy` [documentation](https://scanpy.readthedocs.io/en/latest/api/index.html#recipes).

Depending on the options chosen, `prepare` can take a long time to run (a few minutes for datasets with 10-100k cells, up to an hour or more for datasets with >100k cells). If you want `prepare` to run faster we recommend using the `sparse` option and only computing the embedding for `umap`, using a call like this

```
cellxgene prepare dataset.h5ad --output=dataset-processed.h5ad --embedding=umap --sparse
```

To see all options call

```
cellxgene prepare --help
```


## conda and virtual environments

If you use conda and want to create a conda environment for `cellxgene` you can use the following commands

```
conda create --yes -n cellxgene python=3.7
conda activate cellxgene
pip install cellxgene
```

Or you can create a virtual environment by using

```
ENV_NAME=cellxgene
python3.7 -m venv ${ENV_NAME}
source ${ENV_NAME}/bin/activate
pip install cellxgene
```

## docker

We have included a dockerfile to conveniently run cellxgene from docker.

1. Build the image `docker build . -t cellxgene`
2. Run the container and mount data `docker run -v "$PWD/example-dataset/:/data/" -p 5005:5005 cellxgene launch --host 0.0.0.0 data/pbmc3k.h5ad`
   - You will need to use --host 0.0.0.0 to have the container listen to incoming requests from the browser
