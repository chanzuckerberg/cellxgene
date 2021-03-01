---
title: Install
subtitle: Install
layout: default
---

# Installing cellxgene

Cellxgene has two parts:

- [`cellxgene`](launch) is the main explorer application, which takes an already-processed `h5ad` file as input. This is installed by default.
- [`cellxgene prepare`](prepare) provides auxiliary functionality for preparing your dataset. This is _not_ installed by default.

## Requirements

You'll need **python 3.6+** (Note: At this moment, **python 3.9** is not supported by cellxgene) and an up-to-date version of **Google Chrome**.
The web UI is tested on OSX and Windows using Chrome, and the python CLI is tested on OSX and Ubuntu (via WSL/Windows).
It should work on other platforms, but if you run into trouble let us know.

[Python.org](https://www.python.org/downloads/) has help on installing a recent
version of Python, including the pip package manager. Chrome is available at
[Google.com/chrome](https://google.com/chrome).

## Basic install using pip

To install the `cellxgene` explorer alone, run:

```
pip install cellxgene
```

To install `cellxgene` and the optional `cellxgene prepare`, run:

```
pip install cellxgene[prepare]
```

_Note: if the aforementioned optional `prepare` package installation fails, you can also install these packages directly:_

```
pip install scanpy>=1.3.7 python-igraph louvain>=0.6
```

_On various Linux platforms, you may also need to install build dependencies first:_

```
sudo apt-get install build-essential python-dev
pip install scanpy>=1.3.7 python-igraph louvain>=0.6
```

If you already have `cellxgene` installed, you can update to the most recent version by running:

```
pip install cellxgene --upgrade
```

## Using a conda environment

To install `cellxgene` alone, run:

```
conda create --yes -n cellxgene python=3.7
conda activate cellxgene
pip install cellxgene
```

To install `cellxgene` and the optional `cellxgene prepare`, run:

```
conda create --yes -n cellxgene python=3.7
conda activate cellxgene
pip install cellxgene[prepare]
```

## Using a virtual environment

To install `cellxgene` alone, run:

```
ENV_NAME=cellxgene
python3.7 -m venv ${ENV_NAME}
source ${ENV_NAME}/bin/activate
pip install cellxgene
```

To install `cellxgene` and `cellxgene prepare`, run:

```
ENV_NAME=cellxgene
python3.7 -m venv ${ENV_NAME}
source ${ENV_NAME}/bin/activate
pip install cellxgene[prepare]
```

## Using docker

Build the image

```
docker build . -t cellxgene
```

Run the container and mount data (change data location, `--port` and `--host` parameters as needed)

```
docker run -v "$PWD/example-dataset/:/data/" -p 5005:5005 cellxgene launch --host 0.0.0.0 data/pbmc3k.h5ad
```

You will need to use `--host 0.0.0.0` to have the container listen to incoming requests from the browser
