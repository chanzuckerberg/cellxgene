# Locust Load Test

This directory contains scripts to load test cellxgene's backend. It
primary simulates initial data loading and expression data fetch, which
are the most common data routes. It currently does not include tests
for differential expression or re-clustering routes.

## Prerequisites

You need:

- Python 3.6+, and pip
- cellxgene installed
- install the locust dependencies in `requirements-locust.txt`

## To test

1. Choose to run cellxgene in either single dataset or data root mode.
2. Edit config.py to indicate which datasets to load:
   - in single dataset mode, just set `DataSets=[""]`
   - in dataroot (multi-dataset) mode, add the route names, eg, `DataSets=['foo.cxg', 'bar.cxg']`
3. Launch cellxgene in the appropriate mode
4. launch locust, specifying the correct --host argument
5. point your web browser to the locust http server, usually `http://localhost:8089/`

### Single dataset mode

- Edit config.py and set `DataSets=[""]`
- in a shell, run `cellxgene launch somefile.h5ad`
- launch locust in another shell, `locust --host http://localhost:5005/` (or wherever you are running cellxgene)
- point a browser to the locust port, usually http://localhost:8089/
- run test

### Multi-dataset mode

- Edit config.py and set `DataSets=["datapath1", ...]`
- in a shell, run `cellxgene launch --dataroot path`

The remainder of the steps are same as single dataset.
