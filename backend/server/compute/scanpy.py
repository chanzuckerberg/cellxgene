import importlib
import numpy as np
from anndata import AnnData
"""
Wrapper for various scanpy modules.  Will raise NotImplementedError if the scanpy
module is not installed/available
"""


def get_scanpy_module():
    try:
        sc = importlib.import_module("scanpy")
        # Future: we could enforce versions here, eg, lookat sc.__version__
        return sc
    except ModuleNotFoundError as e:
        raise NotImplementedError("Please install scanpy to enable re-embedding.") from e
    except Exception as e:
        # will capture other ImportError corner cases
        raise NotImplementedError() from e

def get_scanpy_external_module():
    try:
        sc = importlib.import_module("scanpy")
        # Future: we could enforce versions here, eg, lookat sc.__version__
        return sc.external
    except ModuleNotFoundError as e:
        raise NotImplementedError("Please install scanpy to enable re-embedding.") from e
    except Exception as e:
        # will capture other ImportError corner cases
        raise NotImplementedError() from e

def get_samalg_module():
    try:
        samalg = importlib.import_module("samalg")
        # Future: we could enforce versions here, eg, lookat sc.__version__
        return samalg.SAM
    except ModuleNotFoundError as e:
        raise NotImplementedError("Please install sam-alogirthm to enable SAM during re-embedding") from e
    except Exception as e:
        # will capture other ImportError corner cases
        raise NotImplementedError() from e

def scanpy_umap(adata, obs_mask=None, reembedParams = {}, pca_options={}, neighbors_options={}, umap_options={}):
    """
    Given adata and an obs mask, return a new embedding for adata[obs_mask, :]
    as an ndarray of shape (len(obs_mask), N), where N>=2.

    Do NOT mutate adata.
    """

    # backed mode is incompatible with the current implementation
    if adata.isbacked:
        raise NotImplementedError("Backed mode is incompatible with re-embedding")

    # safely get scanpy module, which may not be present.
    sc = get_scanpy_module()

    # https://github.com/theislab/anndata/issues/311
    obs_mask = slice(None) if obs_mask is None else obs_mask
    adata = adata[obs_mask, :].copy()

    for k in list(adata.obsm.keys()):
        del adata.obsm[k]
    for k in list(adata.uns.keys()):
        del adata.uns[k]
    
    doPreprocess = reembedParams.get("doPreprocess",False)
    minCountsCF = reembedParams.get("minCountsCF",0)
    minGenesCF = reembedParams.get("minGenesCF",0)
    minCellsGF = reembedParams.get("minCellsGF",0)
    maxCellsGF = reembedParams.get("maxCellsGF",100)
    minCountsGF = reembedParams.get("minCountsGF",0)
    doSAM = reembedParams.get("doSAM",False)
    nTopGenesHVG = reembedParams.get("nTopGenesHVG",2000)
    nBinsHVG = reembedParams.get("nBins",20)
    doBatch = reembedParams.get("doBatch",False)
    batchMethod = reembedParams.get("batchMethod","Scanorama")
    batchKey = reembedParams.get("batchKey","")
    scanoramaKnn = reembedParams.get("scanoramaKnn",20)
    scanoramaSigma = reembedParams.get("scanoramaSigma",15)
    scanoramaAlpha = reembedParams.get("scanoramaAlpha",0.1)
    scanoramaBatchSize = reembedParams.get("scanoramaBatchSize",5000)
    bbknnNeighborsWithinBatch = reembedParams.get("bbknnNeighborsWithinBatch",3)
    numPCs = reembedParams.get("numPCs",150)
    pcaSolver = reembedParams.get("pcaSolver","randomized")
    neighborsKnn = reembedParams.get("neighborsKnn",20)
    neighborsMethod = reembedParams.get("neighborsMethod","umap")
    distanceMetric = reembedParams.get("distanceMetric","cosine")
    nnaSAM = reembedParams.get("nnaSAM",50)
    weightModeSAM = reembedParams.get("weightModeSAM","dispersion")
    umapMinDist = reembedParams.get("umapMinDist",0.1)
    logTransform = reembedParams.get("logTransform",False)
    scaleData = reembedParams.get("scaleData",False)
    dataLayer = reembedParams.get("dataLayer","X")
    sumNormalizeCells = reembedParams.get("sumNormalizeCells",False)

    if doPreprocess:
        if adata.raw is not None:
            adata_raw = AnnData(X=adata.raw.X)
            adata_raw.var_names = adata.raw.var_names
            adata_raw.obs_names = adata.obs_names
            adata_raw.obs = adata.obs
            for key in adata.var.keys():
                adata_raw.var[key] = adata.var[key]
        else:
            adata_raw = adata
        
        # sc.pp.filter_cells(adata_raw,min_counts=minCountsCF,min_genes=minGenesCF)
        sc.pp.filter_genes(adata_raw, min_counts=minCountsGF)
        sc.pp.filter_genes(adata_raw, min_cells=minCellsGF/100*adata_raw.shape[0])
        sc.pp.filter_genes(adata_raw, max_cells=maxCellsGF/100*adata_raw.shape[0])
        
        if not doSAM:
            batchKey = None if (batchKey == "" or not doBatch) else batchKey
            sc.pp.highly_variable_genes(adata_raw,batch_key=batchKey,flavor='seurat_v3',n_top_genes=nTopGenesHVG, n_bins=nBinsHVG)
    else:
        adata_raw = adata
        
    if doSAM:
        SAM = get_samalg_module()
        sam=SAM(counts = adata_raw, inplace=True)
        sn = "cell_median" if sumNormalizeCells else None
        norm = "log" if logTransform else None
        sam.preprocess_data(sum_norm=sn,thresh_low=0,thresh_high=1.0, norm=norm, min_expression=0,filter_genes=False)
        # do yo SAM
    else:
        # do yo scanpy
        if sumNormalizeCells:
            sc.pp.normalize_total(adata_raw)
        if logTransform:
            sc.pp.log1p(adata_raw)


    
    if not doSAM:
        X = adata_raw.X
        if dataLayer != "X":
            adata_raw.X = adata_raw.layers[dataLayer]
        if scaleData:
            print('Scaling')
            sc.pp.scale(adata_raw,max_value=10)
        sc.pp.pca(adata_raw,n_comps=min(adata.n_vars - 1, numPCs), svd_solver=pcaSolver)
        adata_raw.X = X
    else:
        X = sam.adata.X
        if dataLayer != "X":
            sam.adata.X = sam.adata.layers[dataLayer]        
        preprocessing = "StandardScaler" if scaleData else "Normalizer"
        sam.run(projection=None,weight_mode=weightModeSAM,preprocessing=preprocessing,distance=distanceMetric,num_norm_avg=nnaSAM)
        sam.adata.X = X        
        adata_raw=sam.adata

    if doBatch:
        # do yo batch correction
        sce = get_scanpy_external_module()
        if doSAM:
            adata_batch = sam.adata
        else:
            adata_batch = adata_raw
        
        if batchMethod == "Harmony":
            sce.pp.harmony_integrate(adata_batch,batchKey,adjusted_basis="X_pca")
        elif batchMethod == "BBKNN":
            metric = "angular" if distanceMetric in ["correlation","cosine"] else "euclidean"
            sce.pp.bbknn(adata_batch, batch_key=batchKey, metric=metric, n_pcs=numPCs)
        elif batchMethod == "Scanorama":
            sce.pp.scanorama_integrate(adata_batch, batchKey, basis='X_pca', adjusted_basis='X_pca',
                                    knn=scanoramaKnn, sigma=scanoramaSigma, alpha=scanoramaAlpha,
                                    batch_size=scanoramaBatchSize)
        if doSAM:
            sam.adata = adata_batch
        else:
            adata_raw = adata_batch
    
    if not doSAM or doSAM and batchMethod == "BBKNN":
        if not doBatch or doBatch and batchMethod != "BBKNN":
            sc.pp.neighbors(adata_raw, n_neighbors=neighborsKnn, use_rep="X_pca",method=neighborsMethod, metric=distanceMetric)    
        sc.tl.umap(adata_raw, min_dist=umapMinDist)
    else:
        sam.run_umap(metric=distanceMetric,min_dist=umapMinDist)
        adata_raw.obsm['X_umap'] = sam.adata.obsm['X_umap']
        adata_raw.obsp['connectivities'] = sam.adata.obsp['connectivities']
        
    umap = adata_raw.obsm["X_umap"]
    result = np.full((obs_mask.shape[0], umap.shape[1]), np.NaN)
    result[obs_mask] = umap
    return result, adata_raw.obsp['connectivities']
