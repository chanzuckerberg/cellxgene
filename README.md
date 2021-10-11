<img src="./docs/cellxgene-logo.png" width="300">

# For the Data Science team
This fork implements some of the key features that have been highly requested by the data science team at CZBiohub.

Features include:
- Hotkeys (SHIFT+? to see a tooltip describing all available  hotkeys)
- End-to-end interactive analysis and reembedding, with new embeddings hierarchically organized.
- LIDAR graph interaction mode (the airplane) - Show an interactive tooltip describing the cells underneath your cursor. Very helpful for the color impaired or for large datasets with hundreds of labels.
- Sankey plots
- Leiden clustering
- Label fusion and deletion
- Interactive selection of data layer for expression visualization

### Installation
The easiest way is to install from the precompiled pip distribution located in `dist/`.
Clone this repo, activate your virtual environment of choice (`conda` or `venv`, typically), and install with:

```
pip install excellxgene
```

Ping me on slack (@Alec) if you have any questions!
