import importlib
import numpy as np
import pdb

"""
Wrapper for various scanpy modules.  Will raise NotImplementedError if the scanpy
module is not installed/available
"""


def get_scanpy_module():
    try:
        sc = importlib.import_module("scanpy")
        # Future: we could enforce versions here, eg, lookat sc.__version__
        return sc
    except ImportError:
        print(f"Caught import error lalalala {str(e)}")
        pdb.post_mortem()
        return sc
    except ModuleNotFoundError as e:
        raise NotImplementedError("Please install scanpy to enable UMAP re-embedding") from e
    except Exception as e:
        # will capture other ImportError corner cases
        print(f"Caught exception error dododood {str(e)}")
        pdb.post_mortem()
        raise NotImplementedError() from e


def scanpy_umap(adata, obs_mask=None, pca_options={}, neighbors_options={}, umap_options={}):
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

    sc.pp.pca(adata, zero_center=None, n_comps=min(adata.n_vars - 1, 50), **pca_options)
    sc.pp.neighbors(adata, **neighbors_options)
    sc.tl.umap(adata, **umap_options)

    umap = adata.obsm["X_umap"]
    result = np.full((obs_mask.shape[0], umap.shape[1]), np.NaN)
    result[obs_mask] = umap
    return result
