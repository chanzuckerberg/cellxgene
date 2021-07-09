<img src="./docs/cellxgene-logo.png" width="300">

# For the Data Science team
This fork implements some of the key features that have been highly requested by the data science team at CZBiohub.

Features include:
- Hotkeys (SHIFT+? to see a tooltip describing all available  hotkeys)
- End-to-end interactive analysis and reembedding.

### Installation
The easiest way is to install from the precompiled pip distribution located in `dist/`.
Clone this repo, activate your virtual environment of choice (`conda` or `venv`, typically), and install with:

```
pip install dist/cellxgene-latest.tar.gz sam-algorithm==0.8.6 hnswlib bbknn scanorama
```
Ping me on slack (@Alec) if you have any questions!
