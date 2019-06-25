### How to set up a testing environment for changes related to web hosting.

We often get PRs related to someone using a server to host cellxgene externally or on a local network (ex. https://github.com/chanzuckerberg/cellxgene/pull/568 ). Here is how you can test these changes locally.

We are going to run docker containers for cellxgene and an apache server running a reverse proxy on a local docker network. We run the cellxgene container without exposing any ports so that we cannot access it directly, only through the apache server. We can also update our cellxgene Dockerfile so that we can install a local build instead of having to deploy to pypi.

1 Create Docker network, this allows the containers to communicate with each other.

```
docker network create cxg
```

2 Create and run cellxgene container

(optional) To install cellxgene from the local codebase

a Create sdist file
`make pydist`

b Update Dockerfile to install from dist

```
FROM ubuntu:bionic

ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8
COPY [ "dist/",    "/cellxgene/dist/" ]

RUN apt-get update && \
    apt-get install -y build-essential libxml2-dev python3-dev python3-pip zlib1g-dev && \
    pip3 install /cellxgene/dist/cellxgene-0.5.1.tar.gz

ENTRYPOINT ["cellxgene"]
```

(required) Build container
`docker build . -t cellxgene`

3 Create the proxy container

In a separate directory create these two files

Dockerfile

```
FROM rgoyard/apache-proxy:latest
ADD proxy.conf /conf/
```

proxy.conf

```
ProxyPass "/data/" http://cellxgene:5005/
ProxyPassReverse "/data/" http://cellxgene:5005/
```

Build the container
`docker build -t proxy .`

4 Run containers and attach to network

```
docker run -d -p 80:80 --network cxg --name proxy proxy
docker run -v "$PWD/example-dataset/:/data/" --name cellxgene --network cxg cellxgene launch --host 0.0.0.0 data/pbmc3k.h5ad
```

5 Go to served site

http://localhost/data/
