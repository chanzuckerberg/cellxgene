---
layout: default
title: FAQ
description: Frequently Asked Questions
---

# Data formatting

#### What file formats can I use with _cellxgene_?

Currently, you can go straight into `cellxgene launch` with your own analyzed data in h5ad format, after you have performed dimenstionality reduction (tsne, umap) and clustering (louvain).

If your data is in a different format, and/or you still need to perform dimensionality reduction and clustering, `cellxgene` can do that for you with the `prepare` command. `cellxgene prepare` runs `scanpy` under the hood and can read in any format that is currently supported by `scanpy` (including mtx, loom, and more listed [here](https://scanpy.readthedocs.io/en/latest/api/index.html#reading)).

To add `cellxgene prepare` to your cellxgene installation run `pip install cellxgene[prepare]`. 

The output of `cellxgene prepare` is a h5ad file with your computed clusters and tsne/umap projections that can be used in `cellxgene launch`.

#### I have a directory of 10X-Genomics data with _mtx_ files and I've never used _scanpy_, can I use _cellxgene_?

Yep! This should only take a couple steps. We'll assume your data is in a folder called `data/` and you've successfully installed `cellxgene` with the `prepare` packages as described above. Just run

```
cellxgene prepare data/ --output=data-processed.h5ad --layout=umap
```

Depending on the size of the dataset, this may take some time. Once it's done, call

```
cellxgene launch data-processed.h5ad --layout=umap --open
```

And your web browser should open with an interactive view of your data.

#### I have extra metadata that I want to add to my dataset

Currently this is not supported directly, but you should be able to do this yourself using `scanpy`. For example, this [notebook](https://github.com/falexwolf/fun-analyses/blob/master/tabula_muris/tabula_muris.ipynb) shows adding the contents of a `csv` file with metadata to an `anndata` object. For now, you could do this manually on your data in the same way and then save out the result before loading into `cellxgene`.

#### What part of the _anndata_ objects does cellxgene pull in for visualization?

- `.obs` and `.var` annotations are use to extract metadata for filtering
- `.X` is used to display expression (histograms, scatterplot & colorscale) and to compute differential expression
- `.obsm` is used for layout. If an embedding has more than two components, the first two will be used for visualization.

#### I have a BIG dataset - how can I make cellxgene run as fast as possible?

If your dataset requires gigabytes of disk space, you may need to select an appropriate storage format in order to effectively utilize `cellxgene`. Tips and tricks:

- `cellxgene` is optimized for columnar data access. For large datasets, format the expression matrix (`.X`) as either a [SciPy CSC sparse matrix](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csc_matrix.html) or a dense Numpy array (whichever creates a smaller `h5ad` file). If you are using `cellxgene prepare`, include the `--sparse` flag to ensure `.X` is formatted as a CSC sparse matrix (by default, `.X` will be a dense matrix).
- `cellxgene` start time is directly proportional to `h5ad` file size and the speed of your file system. Expect that large (eg, million cell) datasets will take minutes to load, even on relatively fast computers with a high performance local hard drive. Once loaded, exploring metadata should still be quick.
- If your dataset size exceeds the size of memory (RAM) on the host computer, differential expression calculations will be extremely slow (or fail, if you run out of virtual memory).

# Algorithms

#### How are you computing and sorting differential expression results?

We use a [Welch's _t_-test](https://en.wikipedia.org/wiki/Welch%27s_t-test) implementation including the same variance overestimation correction as used in `scanpy`. We sort the `tscore` to identify the top N genes, and then filter to remove any that fall below a cutoff log fold change value, which can help remove spurious test results. The default threshold is `0.01` and can be changed using the option `--diffexp-lfc-cutoff`.

# Problems, errors, & bugs

#### How do I create a Python environment for _cellxgene_?

If you use conda and want to create a [conda environment](https://conda.io/docs/user-guide/tasks/manage-environments.html) for _cellxgene_ you can use the following commands

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

#### I ran _prepare_ and I'm getting results that look unexpected

You might want to try running one of the preprocessing recipes included with `scanpy` (read more about them [here](https://scanpy.readthedocs.io/en/latest/api/index.html#recipes)). You can specify this with the `--recipe` option, such as

```
cellxgene prepare data/ --output=data-processed.h5ad --recipe=zheng17
```

It should be easy to run `prepare` then call `cellxgene launch` a few times with different settings to explore different behaviors. We may explore adding other preprocessing options in the future.

#### I tried to `pip install cellxgene` and got a weird error I don't understand

This may happen, especially as we work out bugs in our installation process! Please create a new [Github issue](https://github.com/chanzuckerberg/cellxgene/issues), explain what you did, and include all the error messages you saw. It'd also be super helpful if you call `pip freeze` and include the full output alongside your issue.

#### I'm following the developer instructions and get an error about "missing files and directories” when trying to build the client

This is likely because you do not have node and npm installed, we recommend using [nvm](https://github.com/creationix/nvm) if you're new to using these tools.

# Data access

#### Can I use a _s3:_ or _gs:_ URL with `cellxgene launch`?

Yes. Support for S3 and GCS is not enabled by default. If you wish to directly access S3 or GFS, install one or both of the following packages using `pip`:

- [s3fs](https://s3fs.readthedocs.io/en/latest/) for S3 support
- [gcsfs](https://gcsfs.readthedocs.io/en/latest/) for GCS support

For example:

```
pip install s3fs
cellxgene launch s3://mybucket.s3-us-west-2.amazonaws.com/mydata.h5ad
```
