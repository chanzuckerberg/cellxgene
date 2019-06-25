FROM python:3.7

WORKDIR /usr/src/app

RUN pip3 install cellxgene

RUN curl -o pbmc3k.h5ad https://raw.githubusercontent.com/chanzuckerberg/cellxgene/master/example-dataset/pbmc3k.h5ad

expose 5005

CMD ["cellxgene", "launch", "pbmc3k.h5ad", "--host", "0.0.0.0"]
