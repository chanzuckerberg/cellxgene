## dev

<!-- 1. make build-for-server-dev
2. cellxgene launch ./example-dataset/Spatial_Drug.h5ad -->

0. pip install virtualenv
<!-- 1. ./scripts/backend_dev (in $PROJECT_ROOT) -->
1. `cellxgene launch --debug example-dataset/Spatial_Drug.h5ad` (in $PROJECT_ROOT)
2. `make start-frontend` (in $PROJECT_ROOT/client)
3. open `http://localhost:3000`

## build

1. make pydist
2. docker build . -t cellxgene
3. docker run -it -v ./example-dataset:/example-dataset -p 5005:5005 cellxgene launch --host 0.0.0.0 /example-dataset/Spatial_Drug.h5ad

```yaml
services:
  cellxgene:
    container_name: cellxgene
    image: registry.vinoai.net/cellxgene
    command: launch --host 0.0.0.0 /dataset/Spatial_Drug.h5ad
    ports:
      - "5005:5005"
    volumes:
      - ./dataset:/dataset
```
