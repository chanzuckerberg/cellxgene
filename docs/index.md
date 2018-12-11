# cellxgene

cellxgene is an interactive data explorer for single-cell transcriptomics data designed to handle large datasets (1 million cells or more) and integrate with your favorite analysis tools

## getting started

install the package
> `> pip install cellxgene`

preprocess the data for use with cellxgene (optional)
> `> cellxgene --prepare dataset.h5ad -o processed.h5ad`

launch the web app
> `> cellxgene --launch processed.h5ad`

## features



### inspiration and collaboration

We've been heavily inspired by several other related single-cell visualization projects:
* [UCSC Cell Browswer](http://cells.ucsc.edu/)
* [Cytoscape](http://www.cytoscape.org/)
* [Xena](https://xena.ucsc.edu/)
* [ASAP](https://asap.epfl.ch/)
* [Gene Pattern](http://genepattern-notebook.org/)

We were inspired by Mike Bostock and the [crossfilter](https://github.com/crossfilter) team for the design of our filtering implementation.

We have been working closely with the [`scanpy`](https://github.com/theislab/scanpy) team to integrate with their awesome analysis tools. Special thanks to Alex Wolf, Fabian Theis, and the rest of the team for their help during development and for providing an example dataset.

We are eager to explore integrations with other computational backends such as [`Seurat`](https://github.com/satijalab/seurat) or [`Bioconductor`](https://github.com/Bioconductor)

### help and contact

Have questions, suggestions, or comments? You can come hang out with us by joining the [CZI Science Slack](https://join-cziscience-slack.herokuapp.com/) and posting in the `#cellxgene-users` channel. As mentioned above, please submit any feature requests or bugs as [Github issues](https://github.com/chanzuckerberg/cellxgene/issues). We'd love to hear from you!
