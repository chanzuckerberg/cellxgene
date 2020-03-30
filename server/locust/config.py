"""
Locust test config
"""


""" Data routes that will be tested """

# single dataset, for non-dataroot tests
# DataSets = [""]

# multi-dataset, for dataroot tests.  these are varied in size/shape
DataSets = [
    "/d/pbmc3k.cxg",
    "/d/TM_droplet_processed.cxg",
    "/d/pancreas.cxg",
    "/d/immune_bone_marrow_processed.cxg",
    "/d/10X_mouse_13MM_processed.cxg",
    "/d/GSE60361.cxg",
    "/d/Reprogrammed_Dendritic_Cells.cxg",
    "/d/WongAdultRetina.cxg",
]
