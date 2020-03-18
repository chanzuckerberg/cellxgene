"""
Locust test config
"""


""" Data routes that will be tested """

# single dataset, for non-dataroot tests
# DataSets = [""]

# multi-dataset, for dataroot tests.  these are varied in size/shape
DataSets = [
    "/pbmc3k.cxg",
    "/TM_droplet_processed.cxg",
    "/pancreas.cxg",
    "/immune_bone_marrow_processed.cxg",
    "/10x_mouse_13MM_processed.cxg",
]
