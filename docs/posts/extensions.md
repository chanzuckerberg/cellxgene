# cellxgene Extensions

This project was started with the sole goal of empowering the scientific community to explore and understand their data. 
As such, we encourage other scientific tool builders in academia or industry to adopt the patterns, tools, and code from 
this project. All code is freely available for reuse under the [MIT license](https://opensource.org/licenses/MIT).

Before extending cellxgene, we encourage you to reach out to us with ideas or questions. It might be possible that an 
extension could be directly contributed, which would make it available for a wider audience, or that it's on our 
[roadmap](./roadmap.md) and under active development. 

Please note that cellxgene does not have public APIs. Our development may break extensions. We will document changes to the code base but it is advised that extensions pin the version of cellxgene they develop against. 

## Example Reuse & extensions 

#### cellxgene-gateway

[cellxgene Gateway](https://github.com/Novartis/cellxgene-gateway) allows you to use with multiple datasets. It 
displays an index of available h5ad (anndata) files. When a user clicks on a file name, it launches a Cellxgene Server 
instance that loads that particular data file and once it is available proxies requests to that server.

#### cellxgene-VIP (Visualization in Plugin)

[cellxgene-VIP](https://github.com/interactivereport/cellxgene_VIP) enables cellxgene to generate violin, stacked violin, stacked bar, heatmap, volcano, embedding, dot, track, density, 2D density, sankey and dual-gene plot in high-resolution SVG/PNG format. It also performs differential gene expression analysis and provides a Command Line Interface (CLI) for advanced users to perform analysis using python and R.
![cellxgene_VIP](https://interactivereport.github.io/cellxgene_VIP/cellxgene_VIP.png?raw=true "cellxgene_VIP")

#### Galaxy

[Galaxy](https://singlecell.usegalaxy.eu/) is an open source, web-based platform for data intensive biomedical research. cellxgene can be accessed within Galaxy to view analyzed datasets. 
See also the relevant [publication](https://www.biorxiv.org/content/10.1101/2020.06.06.137570v1.full.pdf) 

#### Single Cell Portal

The [Single Cell Portal](https://singlecell.broadinstitute.org/single_cell) is a data hosting and visualization service. cellxgene can be embedded as an additional view to complement the visualizations provided by the. 
[Example](https://singlecell.broadinstitute.org/single_cell/study/SCP807/atlas-of-healthy-and-shiv-infected-non-human-primate-lung-and-ileum-ace2-cells). 

#### FASTGenomics
[FASTGenomics](https://www.fastgenomics.org/) is a collaborative research platform that offers easy-to-use data management and analytics to drive single-cell research forward. cellxgene has been embedded in, and can be used to view datasets processed via . 
See also this 
[Example](https://beta.fastgenomics.org/datasets/detail-dataset-952687f71ef34322a850553c4a24e82e#Cellxgene), 
Note that account creation is required. 
