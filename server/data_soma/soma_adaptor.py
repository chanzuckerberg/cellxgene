import warnings
import tiledbsc

import server.common.compute.diffexp_generic as diffexp_generic
import server.common.compute.estimate_distribution as estimate_distribution
from server.common.errors import PrepareError, DatasetAccessError
from server.common.colors import convert_anndata_category_colors_to_cxg_category_colors
from server.common.fbs.matrix import encode_matrix_fbs

class SomaAdaptor(DataAdaptor):
    def __init__(self, data_locator, app_config=None, dataset_config=None):
        super().__init__(data_locator, app_config, dataset_config)
        self.data = None
        self.X_approximate_distribution = None
        self._load_data(data_locator)
        self._validate_and_initialize()

    @staticmethod
    def pre_load_validation(data_locator):
        # Using same logic as anndata_adaptor
        if data_locator.islocal():
        # if data locator is local, apply file system conventions and other "cheap"
        # validation checks.  If a URI, defer until we actually fetch the data and
        # try to read it.  Many of these tests don't make sense for URIs (eg, extension-
        # based typing).
        if not data_locator.exists():
            raise DatasetAccessError("does not exist")
        if not data_locator.isfile():
            raise DatasetAccessError("is not a file")

    @staticmethod
    def open(data_locator, app_config, dataset_config):
        return SomaAdaptor(data_locator, app_config, dataset_config)

    @staticmethod
    def file_size(data_locator):
        return data_locator.size() if data_locator.islocal() else 0

    def get_name(self):
        """return a string name for this data adaptor"""
        return "cellxgene SOMA adaptor version"

    def get_library_versions(self):
        """return a dictionary of library name to library versions"""
        return dict(tiledbsc=str(tiledbsc.__version__))

    def get_embedding_names(self):
        """
        Return a list of pre-computed embedding names.

        function:
            a) generate list of default layouts
            b) validate layouts are legal.  remove/warn on any that are not
            c) cap total list of layouts at global const MAX_LAYOUTS
        """
        # load default layouts from the data.
        layouts = self.dataset_config.embeddings__names

        if layouts is None or len(layouts) == 0:
            layouts = [key[2:] for key in self.data.obsm.keys() if type(key) == str and key.startswith("X_")]

        # remove invalid layouts
        valid_layouts = []
        obsm_keys = self.data.obsm.keys()
        for layout in layouts:
            layout_name = f"X_{layout}"
            if layout_name not in obsm_keys:
                warnings.warn(f"Ignoring unknown layout name: {layout}.")
            elif not self._is_valid_layout(self.data.obsm[layout_name]):
                warnings.warn(f"Ignoring layout due to malformed shape or data type: {layout}")
            else:
                valid_layouts.append(layout)

        if len(valid_layouts) == 0:
            raise PrepareError("No valid layout data.")

        # cap layouts to MAX_LAYOUTS
        return valid_layouts[0:MAX_LAYOUTS]

    def get_embedding_array(self, ename, dims=2):
        """return an numpy array for the given pre-computed embedding name."""
        full_embedding = self.data.obsm[f"X_{ename}"]
        return full_embedding[:, 0:dims]

    def get_X_array(self, obs_mask=None, var_mask=None):
        """return the X array, possibly filtered by obs_mask or var_mask.
        the return type is either ndarray or scipy.sparse.spmatrix."""
        arr = self.data.X.data.df().to_numpy()
        return arr[obs_mask, var_mask]

    def get_X_approximate_distribution(self) -> XApproximateDistribution:
        """return the approximate distribution of the X matrix."""
        if self.X_approximate_distribution is None:
            """Not yet evaluated."""
            assert self.dataset_config.X_approximate_distribution == "auto"
            m = self.data.X.data.df().to_numpy()
            self.X_approximate_distribution = estimate_distribution.estimate_approximate_distribution(m)

        return self.X_approximate_distribution

    def get_shape(self):
        """return the shape of the data matrix X"""
        return self.data.X.shape()

    def query_var_array(self, term_var):
        return getattr(self.data.var.df(), term_var)

    def query_obs_array(self, term_var):
        return getattr(self.data.obs.df(), term_var)

    def get_colors(self):
        return convert_anndata_category_colors_to_cxg_category_colors(self.data)

    def get_obs_index(self):
        name = self.server_config.single_dataset__obs_names
        if name is None:
            return self.original_obs_index
        else:
            return self.data.obs[name]

    def get_obs_columns(self):
        return self.data.obs.df().columns

    def get_obs_keys(self):
        # return list of keys
        return self.data.obs.keys()

    def get_var_keys(self):
        # return list of keys
        return self.data.var.keys()

    def cleanup(self):
        pass

    @abstractmethod
    def get_schema(self):
        """
        Return current schema
        """
        pass

    def annotation_to_fbs_matrix(self, axis, field=None, uid=None):
        """
        Gets annotation value for each observation
        :param axis: string obs or var
        :param fields: list of keys for annotation to return, returns all annotation values if not set.
        :return: flatbuffer: in fbs/matrix.fbs encoding
        """
        if axis == Axis.OBS:
            if labels is not None and not labels.empty:
                df = self.data.obs.df().join(labels, self.parameters.get("obs_names"))
            else:
                df = self.data.obs.df()
        else:
            df = self.data.var.df()

        if fields is not None and len(fields) > 0:
            df = df[fields]
        return encode_matrix_fbs(df, col_idx=df.columns)

    def compute_diffexp_ttest(self, maskA, maskB, top_n, lfc_cutoff):
        if top_n is None:
            top_n = self.dataset_config.diffexp__top_n
        if lfc_cutoff is None:
            lfc_cutoff = self.dataset_config.diffexp__lfc_cutoff
        return diffexp_generic.diffexp_ttest(self, maskA, maskB, top_n, lfc_cutoff)

    def _load_data(self, data_locator):
        pass

    def _validate_and_initialize(self):
        pass