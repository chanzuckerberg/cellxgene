---
layout: default
title: FAQ
description: Frequently Asked Questions
---


# data formatting

## Someone sent me a directory of `10X-Genomics` data with a `mtx` file and I've never used `scanpy`, can I use `cellxgene`?

Yep! This should only take a couple steps. We'll assume your data is in a folder called `data/` and you've successfully installed `cellxgene` with the `louvain` packages as described above. Just run

```
cellxgene prepare data/ --output=data-processed.h5ad --layout=umap
```

Depending on the size of the dataset, this may take some time. Once it's done, call

```
cellxgene launch data-processed.h5ad --layout=umap --open
```

And your web browser should open with an interactive view of your data.

## In my `prepare` command I received the following error `Warning: louvain module is not installed, no clusters will be calculated. To fix this please install cellxgene with the optional feature louvain enabled`

Louvain clustering requires additional dependencies that are somewhat complex, so we don't include them by default. For now, you need to specify that you want these packages by using

```
pip install cellxgene[louvain]
```

## I ran `prepare` and I'm getting results that look unexpected

You might want to try running one of the preprocessing recipes included with `scanpy` (read more about them [here](https://scanpy.readthedocs.io/en/latest/api/index.html#recipes)). You can specify this with the `--recipe` option, such as

```
cellxgene prepare data/ --output=data-processed.h5ad --recipe=zheng17
```

It should be easy to run `prepare` then call `cellxgene launch` a few times with different settings to explore different behaviors. We may explore adding other preprocessing options in the future.

## I have extra metadata that I want to add to my dataset

Currently this is not supported directly, but you should be able to do this manually using `scanpy`. For example, this [notebook](https://github.com/falexwolf/fun-analyses/blob/master/tabula_muris/tabula_muris.ipynb) shows adding the contents of a `csv` file with metadata to an `anndata` object. For now, you could do this manually on your data in the same way and then save out the result before loading into `cellxgene`.

## What part of the _anndata_ objects does cellxgene pull in for visualization?

- `.obs` and `.var` annotations are use to extract metadata for filtering
- `.X` is used to display expression (histograms, scatterplot & colorscale) and to compute differential expression
- `.obsm` is used for layout

## When I start _cellxgene_, I get an error `Unexpected HTTP response 500, INTERNAL SERVER ERROR -- Out of range float values are not JSON compliant` in the web UI, or `Warning: JSON encoding failure - suggest trying --nan-to-num command line option` in the CLI. What can I do?

At the moment, _cellxgene_ is unable to transmit floating point NaN or Infinity values to the web UI (due to a limitation on data serialization method in use). We expect to resolve this in a future release, but in the meantime, you can work around this issue by starting cellxgene with the `--nan-to-num` command line option, ie, `cellxgene launch data.h5ad --nan-to-num`.

This option will convert all NaNs to zero, and all positive/negative infinities to the min/max of the data element within which the value was found (eg, +Infinity within an `obs` annotation will be converted to the maximum finite value in that annotation). This option will increase startup time, so we recommend only using it when the dataset contains NaN/Infinities.

# installing and building

## I tried to `pip install cellxgene` and got a weird error I don't understand

This may happen, especially as we work out bugs in our installation process! Please create a new [Github issue](https://github.com/chanzuckerberg/cellxgene/issues), explain what you did, and include all the error messages you saw. It'd also be super helpful if you call `pip freeze` and include the full output alongside your issue.

## I'm following the developer instructions and get an error about "missing files and directories” when trying to build the client

This is likely because you do not have node and npm installed, we recommend using [nvm](https://github.com/creationix/nvm) if you're new to using these tools.

# algorithms

## How are you computing and sorting differential expression results?

We use a [Welch's _t_-test](https://en.wikipedia.org/wiki/Welch%27s_t-test) implementation including the same variance overestimation correction as used in `scanpy`. We sort the `tscore` to identify the top N genes, and then filter to remove any that fall below a cutoff log fold change value, which can help remove spurious test results. The default threshold is `0.01` and can be changed using the option `--diffexp-lfc-cutoff`.
