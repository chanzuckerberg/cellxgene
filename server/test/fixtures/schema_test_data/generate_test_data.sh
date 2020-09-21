#!/bin/bash
wget "https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz"
tar xf "pbmc3k_filtered_gene_bc_matrices.tar.gz"

python3 - <<MERGE_GENES
import os
from scipy.io import mmread, mmwrite
import scipy.sparse
import pandas as pd
from server.converters.schema import gene_symbol

mat = mmread("filtered_gene_bc_matrices/hg19/matrix.mtx").todense()
genes = pd.read_csv("filtered_gene_bc_matrices/hg19/genes.tsv", sep='\t', names=["gene_id", "gene_symbol"])

upgraded_genes = gene_symbol.get_upgraded_var_index(pd.DataFrame(index=genes["gene_symbol"]))
df = pd.DataFrame(data=mat, index=upgraded_genes).T
merged = df.sum(axis=1, level=0, skipna=False)

os.makedirs("merged")
merged.columns.to_frame().to_csv("merged/genes.tsv", index=False, header=False)
mmwrite("merged/matrix.mtx", scipy.sparse.coo_matrix(merged).T)
MERGE_GENES

cp "filtered_gene_bc_matrices/hg19/barcodes.tsv" "merged/barcodes.tsv"
awk '{print $1"\t"$1}' merged/genes.tsv > genes_tmp.tsv; mv genes_tmp.tsv merged/genes.tsv

echo -e "\n\n\nRunning tutorial on original\n\n\n"
Rscript - <<TUTORIAL
library(Seurat)

pbmc.data <- Read10X(data.dir = "filtered_gene_bc_matrices/hg19/")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.features = 200)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:10)
saveRDS(pbmc, file = "./seurat_tutorial.rds")
TUTORIAL

echo -e "\n\n\nRunning tutorial on merged\n\n\n"
Rscript - <<TUTORIAL_MERGED
library(Seurat)

pbmc.data <- Read10X(data.dir = "merged/")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.features = 200)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:10)
saveRDS(pbmc, file = "./seurat_tutorial_merged.rds")
TUTORIAL_MERGED

echo -e "\n\n\nRunning SCTransform on original\n\n\n"
Rscript - <<SCTRANSFORM
library(Seurat)
library(sctransform)

pbmc.data <- Read10X(data.dir = "filtered_gene_bc_matrices/hg19/")
pbmc <- CreateSeuratObject(counts = pbmc.data)
pbmc <- PercentageFeatureSet(pbmc, pattern = "^MT-", col.name = "percent.mt")
pbmc <- SCTransform(pbmc, vars.to.regress = "percent.mt", verbose = FALSE)
pbmc <- RunPCA(pbmc, verbose = FALSE)
pbmc <- RunUMAP(pbmc, dims = 1:30, verbose = FALSE)
pbmc <- FindNeighbors(pbmc, dims = 1:30, verbose = FALSE)
pbmc <- FindClusters(pbmc, verbose = FALSE)
saveRDS(pbmc, file = "./sctransform.rds")
SCTRANSFORM

echo -e "\n\n\nRunning SCTransform on merged\n\n\n"
Rscript - <<SCTRANSFORM_MERGED
library(Seurat)
library(sctransform)

pbmc.data <- Read10X(data.dir = "merged/")
pbmc <- CreateSeuratObject(counts = pbmc.data)
pbmc <- PercentageFeatureSet(pbmc, pattern = "^MT-", col.name = "percent.mt")
pbmc <- SCTransform(pbmc, vars.to.regress = "percent.mt", verbose = FALSE)
pbmc <- RunPCA(pbmc, verbose = FALSE)
pbmc <- RunUMAP(pbmc, dims = 1:30, verbose = FALSE)
pbmc <- FindNeighbors(pbmc, dims = 1:30, verbose = FALSE)
pbmc <- FindClusters(pbmc, verbose = FALSE)
saveRDS(pbmc, file = "./sctransform_merged.rds")
SCTRANSFORM_MERGED

echo -e "\n\n\nConverting\n\n\n"
Rscript - <<SCEASY
library(sceasy)
srt <- readRDS("seurat_tutorial.rds")
sceasy::convertFormat(srt,
                      outFile = "seurat_tutorial.h5ad",
                      from = "seurat",
                      to = "anndata",
                      assay = "RNA",
                      main_layer = "data",
                      transfer_layers = c("data", "counts", "scale.data"),
                      drop_single_values = FALSE)

srt <- readRDS("seurat_tutorial_merged.rds")
sceasy::convertFormat(srt,
                      outFile = "seurat_tutorial_merged.h5ad",
                      from = "seurat",
                      to = "anndata",
                      assay = "RNA",
                      main_layer = "data",
                      transfer_layers = c("data", "counts", "scale.data"),
                      drop_single_values = FALSE)

srt <- readRDS("sctransform.rds")
sceasy::convertFormat(srt,
                      outFile = "sctransform.h5ad",
                      from = "seurat",
                      to = "anndata",
                      assay = "SCT",
                      main_layer = "data",
                      transfer_layers = c("data", "counts", "scale.data"),
                      drop_single_values = FALSE)

srt <- readRDS("sctransform_merged.rds")
sceasy::convertFormat(srt,
                      outFile = "sctransform_merged.h5ad",
                      from = "seurat",
                      to = "anndata",
                      assay = "SCT",
                      main_layer = "data",
                      transfer_layers = c("data", "counts", "scale.data"),
                      drop_single_values = FALSE)
SCEASY
