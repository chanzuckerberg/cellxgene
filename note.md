## dev

1. make build-for-server-dev
2. cellxgene launch ./example-dataset/Spatial_Drug.h5ad

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
