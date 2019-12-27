
from server.app.util.matrix_proxy import MatrixProxyView, ArrayProxyView

"""
AnnData/h5py are inconsistent in the API supported by various types of
X matrices.  Sometimes you get a fully ndarray, sometims a Scipy sparse
matrix, sometimes h5py proxies with a subset of our needed interfaces.

This glue code paves over all of that, providing the core set of methods
that the cellxgene ScanPy driver assumes are in existance.  Put another
way, all of the non-portable assumptions are here.
"""


class ArrayProxyView_anndata_h5py(ArrayProxyView):
    """
    override to handle sparse getitem semantics, which differ
    from numpy.
    """
    def toarray(self):
        """ sadly, sparse indexing doesn't drop dimensions like numpy! """
        arr = self.m[self._index[0], self._index[1]]
        if self._vdim == 0:
            arr = arr.transpose()
        return arr.toarray()[0]


class MatrixProxy_anndata_h5py(MatrixProxyView):
    """
    AnnData sparse array stored in H5AD, or proxies for backed data.
    None of these handle indexing very well, so we plop a proxy on top.
    """
    @classmethod
    def __supports__(cls):
        return ("anndata.h5py.h5sparse.SparseDataset",
                "anndata.h5py.h5sparse.backed_csc_matrix",
                "anndata.h5py.h5sparse.backed_csr_matrix",
                "h5py._hl.dataset.Dataset")

    @classmethod
    def create_array(cls, *args, **kwargs):
        return ArrayProxyView_anndata_h5py(*args, **kwargs)
