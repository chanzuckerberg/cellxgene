<img src="./docs/cellxgene-logo.png" width="300">

# Exploratory CellxGene (ExCellxGene)
This fork implements some of the key features that have been highly requested by the data science team at CZBiohub.

Features include:
- Hotkeys (SHIFT+? to see a tooltip describing all available  hotkeys)
- End-to-end interactive analysis and reembedding, with new embeddings hierarchically organized.
- LIDAR graph interaction mode (the airplane) - Show an interactive tooltip describing the cells underneath your cursor. Very helpful for the color impaired or for large datasets with hundreds of labels.
- Sankey plots
- Leiden clustering
- Label fusion and deletion
- Interactive selection of data layer for expression visualization
- Many other quality-of-life improvements.

### Installation

1. Install miniconda if conda not available already:

```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O ~/miniconda.sh
bash ~/miniconda.sh -b -p $HOME/miniconda
```

2. Create and activate a new environment:

```
conda create -n cxg python=3.7
conda activate cxg
```

3. Install excellxgene with pip (the latest version should be `1.1.2`)
```
pip install excellxgene
```

4. Download the git repository to get the example datasets (assumes git is available, if not install it with conda install -c anaconda git)
```
git clone https://github.com/czbiohub/cellxgene
cd cellxgene
```
Datasets are stored in `example-dataset`

5. Launch cellxgene with:
```
cellxgene launch example-dataset
```


This should launch a cellxgene session with all the datasets in example-datasets/ loaded in.

If you're running excellxgene remotely, please launch with:
```
cellxgene launch example-datasets --host 0.0.0.0
```

Ping me on the Biohub slack (@Alec) if you have any questions!
