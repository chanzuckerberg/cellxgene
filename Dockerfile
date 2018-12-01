FROM ubuntu:bionic

ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

EXPOSE 5005:5005

RUN apt-get update
RUN apt-get install -y build-essential libxml2-dev python3-dev python3-pip zlib1g-dev
RUN pip3 install cellxgene

ENTRYPOINT ["cellxgene"]