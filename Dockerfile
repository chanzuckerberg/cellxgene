FROM python:3.7

WORKDIR /usr/src/app

RUN pip3 install cellxgene

expose 5005
