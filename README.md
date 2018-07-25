# cellxgene

### An experiment in interactive, performant exploration of single cell transcriptomics data.

<img align="right" width="350" height="218" src="./example-dataset/cellxgene-demo.gif" pad="50px">
cellxgene is an open-source exploration of how to bring powerful tools from modern web development to visualize and explore large single-cell transcriptomics datasets.
Started in the context of the [Human Cell Atlas Consortium](https://www.humancellatlas.org/), cellxgene hopes to both **enable scientists** to explore their data and to **equip developers** with scalable and reusable patterns and frameworks for visualizing large scientific datasets.

## Features

* **Visualization at scale:** built with [WebGL](https://www.khronos.org/webgl/), [React](https://reactjs.org/) & [Redux](https://redux.js.org/) to handle visualization of at least 1 million cells.

* **Interactive exploration:** select, cross-filter, and compare subsets of your data with performant indexing and data handling.

* **Flexible API:** the cellxgene client-server model is designed to support a range of existing analysis packages for backend computational tasks (eg scanpy), integrated with client-side visualization via a [REST API](https://restfulapi.net/).


## Installation

**Requirements**
- OS: OSX, Windows, Linux
- python 3.6
- npm
- Google Chrome 

**Clone project**  
`git clone https://github.com/chanzuckerberg/cellxgene.git`  

**Install client**  
`cd cellxgene`  
`./bin/build-client`  

**To use with virtual env for python**  
(optional, but recommended)  
`ENV_NAME=cellxgene`  
`python3 -m venv ${ENV_NAME}`  
`source ${ENV_NAME}/bin/activate`  

**Install server**    
`python3 setup.py install`  

**Run (with demo data)**  
`cellxgene`  
*In google chrome, navigate to the viewer via the web address printed in your console.  
E.g.,* `Running on http://0.0.0.0:5005/`


## Contributing
We warmly welcome contributions from the community. Please submit any bug reports and feature requests through github issues. Please submit any direct contributions via a branch + pull request.

## Inspiration and collaboration
We’ve been inspired by several other related efforts in this space, including: [], and hope to explore collaborations where useful.

## Reuse
This project was started with the sole goal of empowering the scientific community to explore and understand their data. As such, we whole-heartedly encourage other scientific tool builders to adopt the patterns, tools, and code from this project, and reach out to us with ideas or questions using Github Issues or Pull Requests. All code is freely available for reuse under the [MIT license](https://opensource.org/licenses/MIT).

*We thank Alex Wolf for the demo dataset.*
