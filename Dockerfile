FROM fastgenomics/scanpy:1.3.1-p36-v2

VOLUME [ "/data" ]

COPY [ "bin/",    "/cellxgene/bin/" ]
COPY [ "server/", "/cellxgene/server/" ]
COPY [ "client/", "/cellxgene/client/" ]
COPY [ "setup.py", "MANIFEST.in", "README.md", "/cellxgene/" ]

WORKDIR /cellxgene

RUN \
    apk --no-cache add --virtual .builddeps bash nodejs nodejs-npm && \
    bash bin/build-client && \
    rm -rf client/node_modules && \
    sed -i 's/\(anndata\|scanpy\|numpy\|pandas\|scipy\).*//' server/requirements.txt && \
    pip3 install --no-cache-dir . && \
    apk del .builddeps

ENTRYPOINT [ "cellxgene", "scanpy", "/data" ]
