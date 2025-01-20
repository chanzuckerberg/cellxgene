FROM python:3.12.8-bookworm

ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

COPY ./dist/cellxgene-1.3.0.tar.gz /tmp/cellxgene-1.3.0.tar.gz

# RUN apt-get update && \
#     apt-get install -y build-essential libxml2-dev python3-dev python3-pip zlib1g-dev && \
#     pip3 install /tmp/cellxgene-1.3.0.tar.gz

RUN pip3 install /tmp/cellxgene-1.3.0.tar.gz

ENTRYPOINT ["cellxgene"]
