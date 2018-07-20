# cellxgene

A React + Redux web application for exploring large scale single cell RNA sequence data.

### Requirements
- OS: OSX, Windows, Linux
- python 3.6
- npm
- Google Chrome


## Installation

#### clone project

    git clone https://github.com/chanzuckerberg/cellxgene.git

#### install client 

    cd cellxgene
    ./bin/build-client

#### To use with virtual env for python (optional, but recommended)

    ENV_NAME=cellxgene
    python3 -m venv ${ENV_NAME}
    source ${ENV_NAME}/bin/activate

#### install server

    python3 setup.py install

#### commandline help

    cellxgene --help
    # For help with the scanpy engine
    cellxgene scanpy --help

#### run (with demo data)

    cellxgene --title PBMC3K scanpy example-dataset/
   
   
*Thanks to Alex Wolf his help with test data*
