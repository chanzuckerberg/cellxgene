import warnings
from datetime import datetime

import anndata
import numpy as np
from packaging import version
import pandas as pd
import scipy as sp
from pandas.core.dtypes.dtypes import CategoricalDtype
from scipy import sparse
from server_timing import Timing as ServerTiming
import time
import os
from glob import glob
import backend.common.compute.diffexp_generic as diffexp_generic
from backend.common.colors import convert_anndata_category_colors_to_cxg_category_colors
from backend.common.constants import Axis, MAX_LAYOUTS
from backend.server.common.corpora import corpora_get_props_from_anndata
from backend.common.errors import PrepareError, DatasetAccessError, FilterError
from backend.common.utils.type_conversion_utils import get_schema_type_hint_of_array
from backend.server.compute.scanpy import get_scanpy_external_module, AnnData, get_samalg_module, get_scanpy_module
from backend.server.data_common.data_adaptor import DataAdaptor
from backend.common.fbs.matrix import encode_matrix_fbs

anndata_version = version.parse(str(anndata.__version__)).release


def anndata_version_is_pre_070():
    major = anndata_version[0]
    minor = anndata_version[1] if len(anndata_version) > 1 else 0
    return major == 0 and minor < 7


class AnndataAdaptor(DataAdaptor):
    def __init__(self, data_locator, app_config=None, dataset_config=None):
        super().__init__(data_locator, app_config, dataset_config)
        self.data = None
        self._load_data(data_locator)
        self._validate_and_initialize()

    def cleanup(self):
        pass

    @staticmethod
    def pre_load_validation(data_locator):
        if data_locator.islocal():
            # if data locator is local, apply file system conventions and other "cheap"
            # validation checks.  If a URI, defer until we actually fetch the data and
            # try to read it.  Many of these tests don't make sense for URIs (eg, extension-
            # based typing).
            if not data_locator.exists():
                raise DatasetAccessError("does not exist")

    @staticmethod
    def file_size(data_locator):
        return data_locator.size() if data_locator.islocal() else 0

    @staticmethod
    def open(data_locator, app_config, dataset_config=None):
        return AnndataAdaptor(data_locator, app_config, dataset_config)

    def get_corpora_props(self):
        return corpora_get_props_from_anndata(self.data)

    def get_name(self):
        return "cellxgene anndata adaptor version"

    def get_library_versions(self):
        return dict(anndata=str(anndata.__version__))

    @staticmethod
    def _create_unique_column_name(df, col_name_prefix):
        """given the columns of a dataframe, and a name prefix, return a column name which
        does not exist in the dataframe, AND which is prefixed by `prefix`

        The approach is to append a numeric suffix, starting at zero and increasing by
        one, until an unused name is found (eg, prefix_0, prefix_1, ...).
        """
        suffix = 0
        while f"{col_name_prefix}{suffix}" in df:
            suffix += 1
        return f"{col_name_prefix}{suffix}"

    def _alias_annotation_names(self):
        """
        The front-end relies on the existance of a unique, human-readable
        index for obs & var (eg, var is typically gene name, obs the cell name).
        The user can specify these via the --obs-names and --var-names config.
        If they are not specified, use the existing index to create them, giving
        the resulting column a unique name (eg, "name").

        In both cases, enforce that the result is unique, and communicate the
        index column name to the front-end via the obs_names and var_names config
        (which is incorporated into the schema).
        """
        self.original_obs_index = self.data.obs.index

        for (ax_name, var_name) in ((Axis.OBS, "obs"), (Axis.VAR, "var")):
            config_name = f"single_dataset__{var_name}_names"
            parameter_name = f"{var_name}_names"
            name = getattr(self.server_config, config_name)
            df_axis = getattr(self.data, str(ax_name))
            if name is None:
                # Default: create unique names from index
                if not df_axis.index.is_unique:
                    raise KeyError(
                        f"Values in {ax_name}.index must be unique. "
                        "Please prepare data to contain unique index values, or specify an "
                        "alternative with --{ax_name}-name."
                    )
                name = self._create_unique_column_name(df_axis.columns, "name_")
                self.parameters[parameter_name] = name
                # reset index to simple range; alias name to point at the
                # previously specified index.
                df_axis.rename_axis(name, inplace=True)
                df_axis.reset_index(inplace=True)
            elif name in df_axis.columns:
                # User has specified alternative column for unique names, and it exists
                if not df_axis[name].is_unique:
                    raise KeyError(
                        f"Values in {ax_name}.{name} must be unique. " "Please prepare data to contain unique values."
                    )
                df_axis.reset_index(drop=True, inplace=True)
                self.parameters[parameter_name] = name
            else:
                # user specified a non-existent column name
                raise KeyError(f"Annotation name {name}, specified in --{ax_name}-name does not exist.")

    def _create_schema(self):
        if self.data.raw is not None:
            layers = [".raw"]+list(self.data.layers.keys())
        else:
            layers = list(self.data.layers.keys())
        
        if "X" not in layers:
            layers = ["X"] + layers

        self.schema = {
            "dataframe": {"nObs": self.cell_count, "nVar": self.gene_count, "type": str(self.data.X.dtype)},
            "annotations": {
                "obs": {"index": self.parameters.get("obs_names"), "columns": []},
                "var": {"index": self.parameters.get("var_names"), "columns": []},
            },
            "layout": {"obs": []},
            "layers": layers
        }
        for ax in Axis:
            curr_axis = getattr(self.data, str(ax))
            for ann in curr_axis:
                ann_schema = {"name": ann, "writable": True}
                ann_schema.update(get_schema_type_hint_of_array(curr_axis[ann]))
                if ann_schema['type']!='categorical':
                    ann_schema['writable']=False
                self.schema["annotations"][ax]["columns"].append(ann_schema)
        for layout in self.get_embedding_names():
            layout_schema = {"name": layout, "type": "float32", "dims": [f"{layout}_0", f"{layout}_1"]}
            self.schema["layout"]["obs"].append(layout_schema)

    def _refresh_layout_schema(self):
        self.schema["layout"]["obs"] = []
        for layout in self.get_embedding_names():
            layout_schema = {"name": layout, "type": "float32", "dims": [f"{layout}_0", f"{layout}_1"]}
            self.schema["layout"]["obs"].append(layout_schema)            

    def get_schema(self):
        return self.schema

    def _load_data(self, data_locator):
        # as of AnnData 0.6.19, backed mode performs initial load fast, but at the
        # cost of significantly slower access to X data.
        try:
            # there is no guarantee data_locator indicates a local file.  The AnnData
            # API will only consume local file objects.  If we get a non-local object,
            # make a copy in tmp, and delete it after we load into memory.
            with data_locator.local_handle() as lh:
                backed = "r" if self.server_config.adaptor__anndata_adaptor__backed else None

                if os.path.isdir(lh) and len(glob(lh+'/*.gz'))==0:
                    filenames = glob(lh+'/*')
                    adatas = []
                    batch = []
                    for file in filenames:
                        if os.path.isdir(file):
                            backed=False
                    
                    for file in filenames:
                        if os.path.isdir(file):
                            sc = get_scanpy_module()
                            adata = sc.read_10x_mtx(file)
                            filt1,_ = sc.pp.filter_cells(adata,min_counts=100, inplace=False)
                            filt2,_ = sc.pp.filter_cells(adata,min_genes=100, inplace=False)
                            filt = np.logical_and(filt1,filt2)
                            adata = adata[filt].copy()
                        elif file.split('.')[-1] =='csv':
                            adata = sc.read_csv(file) 
                            adata.X = sp.sparse.csc_matrix(adata.X)
                        else:
                            adata = anndata.read_h5ad(file, backed=backed)

                        adatas.append(adata)
                        batch.append([file.split('.h5ad')[0].split('/')[-1]]*adata.shape[0])
                    adata = anndata.concat(adatas,join='inner',axis=0)
                    if "orig.ident" not in adata.obs.keys():
                        key = "orig.ident"
                    else:
                        key = f"orig.ident.{str(hex(int(time.time())))[2:]}"
                    adata.obs[key] = pd.Categorical(np.concatenate(batch))
                elif len(glob(lh+'/*.gz'))>0:
                    sc = get_scanpy_module()
                    adata = sc.read_10x_mtx(lh)
                else:
                    adata = anndata.read_h5ad(lh, backed=backed)
                # as of AnnData 0.6.19, backed mode performs initial load fast, but at the
                # cost of significantly slower access to X data.
                adata.obsm["X_root"] = np.zeros((adata.shape[0],2))
                adata.obs_names_make_unique()
                self.data = adata

        except ValueError:
            raise DatasetAccessError(
                "File must be in the .h5ad format. Please read "
                "https://github.com/theislab/scanpy_usage/blob/master/170505_seurat/info_h5ad.md to "
                "learn more about this format. You may be able to convert your file into this format "
                "using `cellxgene prepare`, please run `cellxgene prepare --help` for more "
                "information."
            )
        except MemoryError:
            raise DatasetAccessError("Out of memory - file is too large for available memory.")
        except Exception:
            raise DatasetAccessError(
                "File not found or is inaccessible. File must be an .h5ad object. "
                "Please check your input and try again."
            )

    def _validate_and_initialize(self):
        if anndata_version_is_pre_070():
            warnings.warn(
                "Use of anndata versions older than 0.7 will have serious issues. Please update to at "
                "least anndata 0.7 or later."
            )

        # var and obs column names must be unique
        if not self.data.obs.columns.is_unique or not self.data.var.columns.is_unique:
            raise KeyError("All annotation column names must be unique.")

        self._alias_annotation_names()
        self._validate_data_types()
        self.cell_count = self.data.shape[0]
        self.gene_count = self.data.shape[1]
        self._create_schema()

        # heuristic
        n_values = self.data.shape[0] * self.data.shape[1]
        if (n_values > 1e8 and self.server_config.adaptor__anndata_adaptor__backed is True) or (n_values > 5e8):
            self.parameters.update({"diffexp_may_be_slow": True})

    def _is_valid_layout(self, arr):
        """return True if this layout data is a valid array for front-end presentation:
        * ndarray, dtype float/int/uint
        * with shape (n_obs, >= 2)
        * with all values finite or NaN (no +Inf or -Inf)
        """
        is_valid = type(arr) == np.ndarray and arr.dtype.kind in "fiu"
        is_valid = is_valid and arr.shape[0] == self.data.n_obs and arr.shape[1] >= 2
        is_valid = is_valid and not np.any(np.isinf(arr)) and not np.all(np.isnan(arr))
        return is_valid

    def _validate_data_types(self):
        # The backed API does not support interrogation of the underlying sparsity or sparse matrix type
        # Fake it by asking for a small subarray and testing it.   NOTE: if the user has ignored our
        # anndata <= 0.7 warning, opted for the --backed option, and specified a large, sparse dataset,
        # this "small" indexing request will load the entire X array. This is due to a bug in anndata<=0.7
        # which will load the entire X matrix to fullfill any slicing request if X is sparse.  See
        # user warning in _load_data().
        X0 = self.data.X[0, 0:1]
        if sparse.isspmatrix(X0) and not sparse.isspmatrix_csc(X0):
            warnings.warn(
                "Anndata data matrix is sparse, but not a CSC (columnar) matrix.  "
                "Performance may be improved by using CSC."
            )
        if self.data.X.dtype != "float32":
            warnings.warn(
                f"Anndata data matrix is in {self.data.X.dtype} format not float32. " f"Precision may be truncated."
            )
        for ax in Axis:
            curr_axis = getattr(self.data, str(ax))
            for ann in curr_axis:
                datatype = curr_axis[ann].dtype
                downcast_map = {
                    "int64": "int32",
                    "uint32": "int32",
                    "uint64": "int32",
                    "float64": "float32",
                }
                if datatype in downcast_map:
                    warnings.warn(
                        f"Anndata annotation {ax}:{ann} is in unsupported format: {datatype}. "
                        f"Data will be downcast to {downcast_map[datatype]}."
                    )
                if isinstance(datatype, CategoricalDtype):
                    category_num = len(curr_axis[ann].dtype.categories)
                    if category_num > 500 and category_num > self.dataset_config.presentation__max_categories:
                        warnings.warn(
                            f"{str(ax).title()} annotation '{ann}' has {category_num} categories, this may be "
                            f"cumbersome or slow to display. We recommend setting the "
                            f"--max-category-items option to 500, this will hide categorical "
                            f"annotations with more than 500 categories in the UI"
                        )

    def annotation_to_fbs_matrix(self, axis, fields=None, labels=None):
        if axis == Axis.OBS:
            if labels is not None and not labels.empty:
                df = self.data.obs.copy()
                for c in labels.columns:
                    df[c] = np.array(list(labels[c]))
            else:
                df = self.data.obs
        else:
            df = self.data.var

        if fields is not None and len(fields) > 0:
            df = df[fields]
        return encode_matrix_fbs(df, col_idx=df.columns)

    def get_embedding_names(self):
        """
        Return pre-computed embeddings.

        function:
            a) generate list of default layouts
            b) validate layouts are legal.  remove/warn on any that are not
            c) cap total list of layouts at global const MAX_LAYOUTS
        """
        # load default layouts from the data.
        layouts = self.dataset_config.embeddings__names

        if layouts is None or len(layouts) == 0:
            layouts = [key[2:] for key in self.data.obsm_keys() if type(key) == str and key.startswith("X_")]

        # remove invalid layouts
        valid_layouts = []
        obsm_keys = self.data.obsm_keys()
        for layout in layouts:
            layout_name = f"X_{layout}"
            if layout_name not in obsm_keys:
                warnings.warn(f"Ignoring unknown layout name: {layout}.")
            elif not self._is_valid_layout(self.data.obsm[layout_name]):
                warnings.warn(f"Ignoring layout due to malformed shape or data type: {layout}")
            else:
                valid_layouts.append(layout)
        # cap layouts to MAX_LAYOUTS
        return valid_layouts[0:MAX_LAYOUTS]

    def get_embedding_array(self, ename, dims=2):
        full_embedding = self.data.obsm[f"X_{ename}"]
        return full_embedding[:, 0:dims]
    
    def compute_leiden(self,name,cName,resolution):
        nnm = self.data.uns.get('N_'+name,None)
        if nnm is None:
            nnm = self.data.obsp['connectivities']
            mask = np.array([True]*nnm.shape[0])
        else:
            mask = self.data.uns['N_'+name+'_mask']
            nnm = nnm[mask][:,mask]        
        X = nnm

        import igraph as ig
        import leidenalg

        adjacency = X
        sources, targets = adjacency.nonzero()
        weights = adjacency[sources, targets]
        if isinstance(weights, np.matrix):
            weights = weights.A1
        g = ig.Graph(directed=True)
        g.add_vertices(adjacency.shape[0])
        g.add_edges(list(zip(sources, targets)))
        try:
            g.es["weight"] = weights
        except BaseException:
            pass

        cl = leidenalg.find_partition(
            g, leidenalg.RBConfigurationVertexPartition, resolution_parameter=resolution,seed=0
        )
        result = np.array(cl.membership)
        clusters = np.array(["unassigned"]*mask.size,dtype='object')
        clusters[mask] = result.astype('str')
        
        self.data.obs[cName] = pd.Categorical(clusters)        
        return result      

    def compute_sankey_df(self, labels, name):
        nnm = self.data.uns.get('N_'+name,None)
        if nnm is None:
            nnm = self.data.obsp['connectivities']
        else:
            mask = self.data.uns['N_'+name+'_mask']
            nnm = nnm[mask][:,mask]
        nnm = (nnm+nnm.T)/2
        
        cl=[]
        clu = []
        for i,c in enumerate(labels):
            cl.append(np.array(['A'+str(i)+'_'+str(x).replace(' ','_').replace('(','_').replace(')','_') for x in c]))
            clu.append(np.unique(cl[-1]))

        ps = []
        cs = []
        kmean = nnm.sum(1).A.mean()
        for i,cl1 in enumerate(cl):
            for cl2 in cl[(i+1):]:
                clu1 = np.unique(cl1)
                clu2 = np.unique(cl2)

                for i, c1 in enumerate(clu1):
                    for j, c2 in enumerate(clu2):
                        if c1.split('_')[1] !='unassigned' and c2.split('_')[1]!='unassigned':
                            val = max(nnm[cl1 == c1,:][:, cl2 == c2].sum(1).A.mean(),
                                nnm[cl2 == c2,:][:, cl1 == c1].sum(1).A.mean())
                            val /= kmean
                            if val > 0.2:
                                cs.append(val)
                                ps.append([c1,c2])
                            
        ps = np.vstack(ps)
        cs = np.array(cs)
        return ps,cs
    
    
    def compute_preprocess(self, reembedParams):
        if 'orig.exprs' not in self.data.layers.keys(): #this is the first time preprocessing :D
            self.data.layers['orig.exprs'] = self.data.X
            
        self.data.X = self.data.layers['orig.exprs'].copy()
        
        adata = self.data

        if adata.isbacked:
            raise NotImplementedError("Backed mode is incompatible with preprocessing.")

        # safely get scanpy module, which may not be present.
        sc = get_scanpy_module()

        cn = np.array(list(adata.obs["name_0"]))
        uns = adata.uns.copy()
        for k in list(adata.uns.keys()):
            del adata.uns[k]
        
        doBatchPrep = reembedParams.get("doBatchPrep",False)
        batchPrepParams = reembedParams.get("batchPrepParams",{})
        batchPrepKey = reembedParams.get("batchPrepKey","")
        batchPrepLabel = reembedParams.get("batchPrepLabel","")

        doPreprocess = reembedParams.get("doPreprocess",False)
        minCountsCF = reembedParams.get("minCountsCF",0)
        minGenesCF = reembedParams.get("minGenesCF",0)
        minCellsGF = reembedParams.get("minCellsGF",0)
        maxCellsGF = reembedParams.get("maxCellsGF",100)
        minCountsGF = reembedParams.get("minCountsGF",0)
        logTransform = reembedParams.get("logTransform",False)
        dataLayer = reembedParams.get("dataLayer","X")
        sumNormalizeCells = reembedParams.get("sumNormalizeCells",False)

            
        filt = np.array([True]*adata.shape[0])

        if doBatchPrep and batchPrepKey != "" and batchPrepLabel != "":
            cl = np.array(list(adata.obs[batchPrepKey]))
            batches = np.unique(cl)
            adatas = []
            cns = []
            for k in batches:
                params = batchPrepParams[batchPrepKey].get(k,{})

                doPreprocess = params.get("doPreprocess",False)
                minCountsCF = params.get("minCountsCF",0)
                minGenesCF = params.get("minGenesCF",0)
                minCellsGF = params.get("minCellsGF",0)
                maxCellsGF = params.get("maxCellsGF",100)
                minCountsGF = params.get("minCountsGF",0)
                logTransform = params.get("logTransform",False)
                dataLayer = params.get("dataLayer","X")
                sumNormalizeCells = params.get("sumNormalizeCells",False)
                
                adata_sub = adata[cl==k].copy()
                adata_sub.obs_names = adata_sub.obs["name_0"]
                if dataLayer == ".raw" and adata_sub.raw is not None:
                    adata_sub_raw = AnnData(X=adata_sub.raw.X)
                    adata_sub_raw.var_names = adata_sub.raw.var_names
                    adata_sub_raw.obs_names = adata_sub.obs_names
                    adata_sub_raw.obs = adata_sub.obs
                    for key in adata_sub.var.keys():
                        adata_sub_raw.var[key] = adata_sub.var[key]
                elif dataLayer == ".raw":
                    adata_sub_raw = adata_sub
                elif dataLayer == "X":
                    adata_sub_raw = adata_sub
                    if dataLayer == "X" and "X" not in adata_sub_raw.layers.keys():
                        adata_sub_raw.layers["X"] = adata_sub_raw.X    
                    adata_sub_raw.X = adata_sub_raw.layers[dataLayer]        
                else:
                    adata_sub_raw = AnnData(X=adata_sub.layers[dataLayer])
                    adata_sub_raw.var_names = adata_sub.raw.var_names
                    adata_sub_raw.obs_names = adata_sub.obs_names
                    adata_sub_raw.obs = adata_sub.obs
                    for key in adata_sub.var.keys():
                        adata_sub_raw.var[key] = adata_sub.var[key]   
                if doPreprocess:
                    filt1,_ = sc.pp.filter_cells(adata_sub_raw,min_counts=minCountsCF, inplace=False)
                    filt2,_ = sc.pp.filter_cells(adata_sub_raw,min_genes=minGenesCF, inplace=False)
                    filt = np.logical_and(filt1,filt2)
                    cns.extend(np.array(list(adata_sub_raw.obs["name_0"]))[filt])
                    target_sum = np.median(np.array(adata_sub_raw.X[filt].sum(1)).flatten())
                    a1,_=sc.pp.filter_genes(adata_sub_raw, min_counts=minCountsGF,inplace=False)
                    a2,_=sc.pp.filter_genes(adata_sub_raw, min_cells=minCellsGF/100*adata_sub_raw.shape[0],inplace=False)
                    a3,_=sc.pp.filter_genes(adata_sub_raw, max_cells=maxCellsGF/100*adata_sub_raw.shape[0],inplace=False)
                    a = a1*a2*a3
                    if sp.sparse.issparse(adata_sub_raw.X):
                        adata_sub_raw.X = adata_sub_raw.X.multiply(a.flatten()[None,:]).tocsc()
                    else:
                        adata_sub_raw.X = adata_sub_raw.X * (a.flatten()[None,:])
                

                    if sumNormalizeCells:
                        sc.pp.normalize_total(adata_sub_raw,target_sum=target_sum)
                    if logTransform:
                        sc.pp.log1p(adata_sub_raw)  
                else: 
                    cns.extend(np.array(list(adata_sub_raw.obs["name_0"])))

                adatas.append(adata_sub_raw)
            adata_raw = anndata.concat(adatas,axis=0,join="inner")
            filt = np.in1d(np.array(list(cn)),np.array(cns))
            temp = adata_raw.obs_names.copy()
            adata_raw.obs_names = adata_raw.obs["name_0"]
            adata_raw = adata_raw[cn]
            adata_raw.obs_names = temp
        else:
            if dataLayer == ".raw" and adata.raw is not None:
                adata_raw = AnnData(X=adata.raw.X)
                adata_raw.var_names = adata.raw.var_names
                adata_raw.obs_names = adata.obs_names
                adata_raw.obs = adata.obs
                for key in adata.var.keys():
                    adata_raw.var[key] = adata.var[key]
            elif dataLayer == ".raw":
                adata_raw = adata.copy()
            elif dataLayer == "X":
                adata_raw = adata.copy()
                if dataLayer == "X" and "X" not in adata_raw.layers.keys():
                    adata_raw.layers["X"] = adata_raw.X    
                adata_raw.X = adata_raw.layers[dataLayer]        
            else:
                adata_raw = AnnData(X=adata.layers[dataLayer])
                adata_raw.var_names = adata.raw.var_names
                adata_raw.obs_names = adata.obs_names
                adata_raw.obs = adata.obs
                for key in adata.var.keys():
                    adata_raw.var[key] = adata.var[key]                
            
            if doPreprocess:
                filt1,_ = sc.pp.filter_cells(adata_raw,min_counts=minCountsCF, inplace=False)
                filt2,_ = sc.pp.filter_cells(adata_raw,min_genes=minGenesCF, inplace=False)
                filt = np.logical_and(filt1,filt2)
                target_sum = np.median(np.array(adata_raw.X[filt].sum(1)).flatten())
                a1,_=sc.pp.filter_genes(adata_raw, min_counts=minCountsGF,inplace=False)
                a2,_=sc.pp.filter_genes(adata_raw, min_cells=minCellsGF/100*adata_raw.shape[0],inplace=False)
                a3,_=sc.pp.filter_genes(adata_raw, max_cells=maxCellsGF/100*adata_raw.shape[0],inplace=False)
                a = a1*a2*a3
                if sp.sparse.issparse(adata_raw.X):
                    adata_raw.X = adata_raw.X.multiply(a.flatten()[None,:]).tocsc()
                else:
                    adata_raw.X = adata_raw.X * (a.flatten()[None,:])
            
                if sumNormalizeCells:
                    sc.pp.normalize_total(adata_raw,target_sum=target_sum)
                if logTransform:
                    sc.pp.log1p(adata_raw) 


        self.data.X = adata_raw.X
        self.data.layers['X'] = adata_raw.X
        self.data.uns=uns
        
        for k in self.data.obsm.keys():
            umap = self.data.obsm[k]
            result = np.full((filt.size, umap.shape[1]), np.NaN)
            result[filt] = umap[filt]
            self.data.obsm[k] = result
        self._create_schema()
        return self.get_schema()

    def compute_embedding(self, method, obsFilter, reembedParams, parentName, embName):
        if Axis.VAR in obsFilter:
            raise FilterError("Observation filters may not contain variable conditions")
        if method != "umap":
            raise NotImplementedError(f"re-embedding method {method} is not available.")
        try:
            shape = self.get_shape()
            obs_mask = self._axis_filter_to_mask(Axis.OBS, obsFilter["obs"], shape[0])
        except (KeyError, IndexError):
            raise FilterError("Error parsing filter")
        with ServerTiming.time("layout.compute"):            
            adata = self.data

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
            scaleData = reembedParams.get("scaleData",False)

 
            if not doSAM:
                try:
                    sc.pp.highly_variable_genes(adata,flavor='seurat_v3',n_top_genes=min(nTopGenesHVG,adata.shape[1]), n_bins=nBinsHVG)                
                    adata = adata[:,adata.var['highly_variable']]                
                except:
                    print('Error during HVG selection - some of your expressions are probably negative.')
                X = adata.X
                if scaleData:
                    sc.pp.scale(adata,max_value=10)
                sc.pp.pca(adata,n_comps=min(adata.n_vars - 1, numPCs), svd_solver=pcaSolver)
                adata.X = X
            else:
                SAM = get_samalg_module()
                sam=SAM(counts = adata, inplace=True)
                X = sam.adata.X
                preprocessing = "StandardScaler" if scaleData else "Normalizer"
                sam.run(projection=None,weight_mode=weightModeSAM,preprocessing=preprocessing,distance=distanceMetric,num_norm_avg=nnaSAM)
                sam.adata.X = X        
                adata=sam.adata

            if doBatch:
                sce = get_scanpy_external_module()
                if doSAM:
                    adata_batch = sam.adata
                else:
                    adata_batch = adata
                
                if batchMethod == "Harmony":
                    sce.pp.harmony_integrate(adata_batch,batchKey,adjusted_basis="X_pca")
                elif batchMethod == "BBKNN":
                    sce.pp.bbknn(adata_batch, batch_key=batchKey, metric=distanceMetric, n_pcs=numPCs, neighbors_within_batch=bbknnNeighborsWithinBatch)
                elif batchMethod == "Scanorama":
                    sce.pp.scanorama_integrate(adata_batch, batchKey, basis='X_pca', adjusted_basis='X_pca',
                                            knn=scanoramaKnn, sigma=scanoramaSigma, alpha=scanoramaAlpha,
                                            batch_size=scanoramaBatchSize)
                if doSAM:
                    sam.adata = adata_batch
                else:
                    adata = adata_batch

            if not doSAM or doSAM and batchMethod == "BBKNN":
                if not doBatch or doBatch and batchMethod != "BBKNN":
                    sc.pp.neighbors(adata, n_neighbors=neighborsKnn, use_rep="X_pca",method=neighborsMethod, metric=distanceMetric)    
                sc.tl.umap(adata, min_dist=umapMinDist)
            else:
                sam.run_umap(metric=distanceMetric,min_dist=umapMinDist)
                adata.obsm['X_umap'] = sam.adata.obsm['X_umap']
                adata.obsp['connectivities'] = sam.adata.obsp['connectivities']
                
            umap = adata.obsm["X_umap"]
            result = np.full((obs_mask.shape[0], umap.shape[1]), np.NaN)
            result[obs_mask] = umap
            X_umap,nnm = result, adata.obsp['connectivities']            

        # Server picks reemedding name, which must not collide with any other
        # embedding name generated by this backend.
        if embName == "":
            embName = f"{method}_{str(hex(int(time.time())))[2:]}"

        if parentName != "":
            parentName+=";;"
        
        name = f"{parentName}{embName}"
        if "X_"+name in self.data.obsm.keys():
            name = f"{name}_{str(hex(int(time.time())))[2:]}"
            
        dims = [f"{name}_0", f"{name}_1"]
        layout_schema = {"name": name, "type": "float32", "dims": dims}
        self.schema["layout"]["obs"].append(layout_schema)

        IXer = pd.Series(index =np.arange(nnm.shape[0]), data = np.where(obs_mask.flatten())[0])
        x,y = nnm.nonzero()
        d = nnm.data
        nnm = sp.sparse.coo_matrix((d,(IXer[x].values,IXer[y].values)),shape=(self.data.shape[0],)*2).tocsr()

        self.data.obsm[f"X_{name}"] = X_umap
        self.data.obsp["connectivities"] = nnm
        self.data.uns[f"N_{name}"] = nnm
        self.data.uns[f"N_{name}_mask"] = np.array(list(obs_mask)).flatten()

        return layout_schema

    def compute_diffexp_ttest(self, maskA, maskB, top_n=None, lfc_cutoff=None):
        if top_n is None:
            top_n = self.dataset_config.diffexp__top_n
        if lfc_cutoff is None:
            lfc_cutoff = self.dataset_config.diffexp__lfc_cutoff
        return diffexp_generic.diffexp_ttest(self, maskA, maskB, top_n, lfc_cutoff)

    def get_colors(self):
        return convert_anndata_category_colors_to_cxg_category_colors(self.data)

    def get_X_array(self, obs_mask=None, var_mask=None):
        if obs_mask is None:
            obs_mask = slice(None)
        if var_mask is None:
            var_mask = slice(None)
        X = self.data.X[obs_mask, var_mask]
        return X

    def get_shape(self):
        return self.data.shape

    def query_var_array(self, term_name):
        return getattr(self.data.var, term_name)

    def query_obs_array(self, term_name):
        return getattr(self.data.obs, term_name)

    def get_obs_index(self):
        name = self.server_config.single_dataset__obs_names
        if name is None:
            return self.original_obs_index
        else:
            return self.data.obs[name]

    def get_obs_columns(self):
        return self.data.obs.columns

    def get_obs_keys(self):
        # return list of keys
        return self.data.obs.keys().to_list()

    def get_var_keys(self):
        # return list of keys
        return self.data.var.keys().to_list()
