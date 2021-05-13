FROM ubuntu:bionic

ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

RUN apt-get update && \
    apt-get install -y build-essential libxml2-dev python3-dev python3-pip zlib1g-dev python3-requests python3-aiohttp && \
    python3 -m pip install --upgrade pip && \
    pip3 install cellxgene

# Temporary workaround for dependency issue. See:
# https://github.com/chanzuckerberg/cellxgene/issues/2197
RUN pip3 install -U flask==1.1.2 flatbuffers==1.12

ENTRYPOINT ["cellxgene"]
