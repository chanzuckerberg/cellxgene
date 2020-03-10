from abc import ABCMeta, abstractmethod
from server_timing import Timing as ServerTiming
import numpy as np
import pandas as pd
from os.path import basename, splitext

from server.data_common.fbs.matrix import encode_matrix_fbs
from server.common.constants import Axis, DEFAULT_TOP_N
from server.common.errors import FilterError, JSONEncodingValueError
from server.compute.diffexp import diffexp_ttest
from server.common.utils import jsonify_numpy
from server.common.app_config import AppFeature, AppConfig


class DataAdaptor(metaclass=ABCMeta):
    """Base class for loading and accessing matrix data"""

    def __init__(self, config):
        # config will normally be a type that inherits from AppConfig.
        # the following is for backwards compatability with tests
        if config is None:
            config = AppConfig()
        elif type(config) == dict:
            config = AppConfig(**config)

        # config is the application configuration
        self.config = config

        # parameters set by this data adaptor based on the data.
        self.parameters = {}

    @staticmethod
    @abstractmethod
    def pre_load_validation(data_locator):
        pass

    @staticmethod
    @abstractmethod
    def open(data_locator, config):
        pass

    @staticmethod
    @abstractmethod
    def file_size(data_locator):
        pass

    @abstractmethod
    def get_name(self):
        """return a string name for this data adaptor"""
        pass

    @abstractmethod
    def get_library_versions(self):
        """return a dictionary of library name to library versions"""
        pass

    @abstractmethod
    def get_embedding_names(self):
        """return a list of pre-computed embedding names"""
        pass

    @abstractmethod
    def get_embedding_array(self, ename, dims=2):
        """return an numpy array for the given pre-computed embedding name."""
        pass

    @abstractmethod
    def compute_embedding(self, method, filter):
        """compute a new embedding on the specified obs subset, and return a
           tuple of (schema, fbs)."""
        pass

    @abstractmethod
    def get_X_array(self, obs_mask=None, var_mask=None):
        """return the X array, possibly filtered by obs_mask or var_mask.
        the return type is either ndarray or scipy.sparse.spmatrix."""
        pass

    @abstractmethod
    def get_shape(self):
        pass

    @abstractmethod
    def query_var_array(self, term_var):
        pass

    @abstractmethod
    def query_obs_array(self, term_var):
        pass

    @abstractmethod
    def get_obs_index(self):
        pass

    @abstractmethod
    def get_obs_columns(self):
        pass

    @abstractmethod
    def cleanup(self):
        pass

    @abstractmethod
    def get_location(self):
        pass

    @abstractmethod
    def get_data_locator(self):
        pass

    def get_about(self):
        return None

    def get_title(self):
        # default to file name
        location = self.get_location()
        if location.endswith("/"):
            location = location[:-1]
        return splitext(basename(location))[0]

    @abstractmethod
    def get_schema(self):
        """
        Return current schema
        """
        pass

    @abstractmethod
    def annotation_to_fbs_matrix(self, axis, field=None, uid=None):
        """
        Gets annotation value for each observation
        :param axis: string obs or var
        :param fields: list of keys for annotation to return, returns all annotation values if not set.
        :return: flatbuffer: in fbs/matrix.fbs encoding
        """
        pass

    def get_features(self, annotations=None):
        """Return list of features, to return as part of the config route"""
        features = [
            AppFeature("/cluster/", method="POST", available=False),
            AppFeature("/layout/obs", method="GET", available=self.get_embedding_names() is not None),
            AppFeature("/layout/obs", method="PUT", available=self.config.enable_reembedding),
            AppFeature("/diffexp/", method="POST", available=not self.config.disable_diffexp),
            AppFeature("/annotations/obs", method="PUT", available=annotations is not None),
        ]
        return features

    def update_parameters(self, parameters):
        parameters.update(self.parameters)

    def _index_filter_to_mask(self, filter, count):
        mask = np.zeros((count,), dtype=np.bool)
        for i in filter:
            if type(i) == list:
                mask[i[0] : i[1]] = True
            else:
                mask[i] = True
        return mask

    def _axis_filter_to_mask(self, axis, filter, count):
        mask = np.ones((count,), dtype=np.bool)
        if "index" in filter:
            mask = np.logical_and(mask, self._index_filter_to_mask(filter["index"], count))
        if "annotation_value" in filter:
            mask = np.logical_and(mask, self._annotation_filter_to_mask(axis, filter["annotation_value"], count))

        return mask

    def _annotation_filter_to_mask(self, axis, filter, count):
        mask = np.ones((count,), dtype=np.bool)
        for v in filter:
            name = v["name"]
            if axis == Axis.VAR:
                anno_data = self.query_var_array(name)
            elif axis == Axis.OBS:
                anno_data = self.query_obs_array(name)

            if anno_data.dtype.name in ["boolean", "category", "object"]:
                values = v.get("values", [])
                key_idx = np.in1d(anno_data, values)
                mask = np.logical_and(mask, key_idx)

            else:
                min_ = v.get("min", None)
                max_ = v.get("max", None)
                if min_ is not None:
                    key_idx = (anno_data >= min_).ravel()
                    mask = np.logical_and(mask, key_idx)
                if max_ is not None:
                    key_idx = (anno_data <= max_).ravel()
                    mask = np.logical_and(mask, key_idx)

        return mask

    def _filter_to_mask(self, filter):
        """
        Return the filter as a row and column selection list.
        No filter on a dimension means 'all'
        """
        shape = self.get_shape()
        var_selector = None
        obs_selector = None
        if filter is not None:
            if Axis.OBS in filter:
                obs_selector = self._axis_filter_to_mask(Axis.OBS, filter["obs"], shape[0])

            if Axis.VAR in filter:
                var_selector = self._axis_filter_to_mask(Axis.VAR, filter["var"], shape[1])

        return (obs_selector, var_selector)

    def check_new_labels(self, labels_df):
        """Check the new annotations labels, then set the labels_df index"""
        if labels_df is None or labels_df.empty:
            return

        labels_df.index = self.get_obs_index()
        if labels_df.index.name is None:
            labels_df.index.name = "index"

        # all labels must have a name, which must be unique and not used in obs column names
        if not labels_df.columns.is_unique:
            raise KeyError(f"All column names specified in user annotations must be unique.")

        # the label index must be unique, and must have same values the anndata obs index
        if not labels_df.index.is_unique:
            raise KeyError(f"All row index values specified in user annotations must be unique.")

        obs_columns = self.get_obs_columns()

        duplicate_columns = list(set(labels_df.columns) & set(obs_columns))
        if len(duplicate_columns) > 0:
            raise KeyError(
                f"Labels file may not contain column names which overlap " f"with h5ad obs columns {duplicate_columns}"
            )

        # labels must have same count as obs annotations
        shape = self.get_shape()
        if labels_df.shape[0] != shape[0]:
            raise ValueError("Labels file must have same number of rows as data file.")

    def data_frame_to_fbs_matrix(self, filter, axis):
        """
        Retrieves data 'X' and returns in a flatbuffer Matrix.
        :param filter: filter: dictionary with filter params
        :param axis: string obs or var
        :return: flatbuffer Matrix

        Caveats:
        * currently only supports access on VAR axis
        * currently only supports filtering on VAR axis
        """

        if axis != Axis.VAR:
            raise ValueError("Only VAR dimension access is supported")

        try:
            obs_selector, var_selector = self._filter_to_mask(filter)
        except (KeyError, IndexError, TypeError, AttributeError) as e:
            raise FilterError(f"Error parsing filter: {e}") from e

        if obs_selector is not None:
            raise FilterError("filtering on obs unsupported")

        X = self.get_X_array(obs_selector, var_selector)
        col_idx = np.nonzero([] if var_selector is None else var_selector)[0]
        return encode_matrix_fbs(X, col_idx=col_idx, row_idx=None)

    def diffexp_topN(self, obsFilterA, obsFilterB, top_n=None):
        """
        Computes the top N differentially expressed variables between two observation sets. If mode
        is "TOP_N", then stats for the top N
        dataframes
        contain a subset of variables, then statistics for all variables will be returned, otherwise
        only the top N vars will be returned.
        :param obsFilterA: filter: dictionary with filter params for first set of observations
        :param obsFilterB: filter: dictionary with filter params for second set of observations
        :param top_n: Limit results to top N (Top var mode only)
        :return: top N genes and corresponding stats
        """
        if Axis.VAR in obsFilterA or Axis.VAR in obsFilterB:
            raise FilterError("Observation filters may not contain variable conditions")
        try:
            shape = self.get_shape()
            obs_mask_A = self._axis_filter_to_mask(Axis.OBS, obsFilterA["obs"], shape[0])
            obs_mask_B = self._axis_filter_to_mask(Axis.OBS, obsFilterB["obs"], shape[0])
        except (KeyError, IndexError) as e:
            raise FilterError(f"Error parsing filter: {e}") from e
        if top_n is None:
            top_n = DEFAULT_TOP_N

        result = diffexp_ttest(self, obs_mask_A, obs_mask_B, top_n, self.config.diffexp_lfc_cutoff)

        try:
            return jsonify_numpy(result)
        except ValueError:
            raise JSONEncodingValueError("Error encoding differential expression to JSON")

    @staticmethod
    def normalize_embedding(embedding):
        """Normalize embedding layout to meet client assumptions.
           Embedding is an ndarray, shape (n_obs, n)., where n is normally 2
        """

        # scale isotropically
        min = embedding.min(axis=0)
        max = embedding.max(axis=0)
        scale = np.amax(max - min)
        normalized_layout = (embedding - min) / scale

        # translate to center on both axis
        translate = 0.5 - ((max - min) / scale / 2)
        normalized_layout = normalized_layout + translate

        normalized_layout = normalized_layout.astype(dtype=np.float32)
        return normalized_layout

    def layout_to_fbs_matrix(self):
        """ same as layout, except returns a flatbuffer """
        """
        return all embeddings as a flatbuffer, using the cellxgene matrix fbs encoding.

        * returns only first two dimensions, with name {ename}_0 and {ename}_1,
          where {ename} is the embedding name.
        * client assumes each will be individually centered & scaled (isotropically)
          to a [0, 1] range.
        * does not support filtering

        """

        embeddings = self.get_embedding_names()
        layout_data = []
        with ServerTiming.time(f"layout.query"):
            for ename in embeddings:
                embedding = self.get_embedding_array(ename, 2)
                normalized_layout = DataAdaptor.normalize_embedding(embedding)
                layout_data.append(pd.DataFrame(normalized_layout, columns=[f"{ename}_0", f"{ename}_1"]))

        with ServerTiming.time(f"layout.encode"):
            if layout_data:
                df = pd.concat(layout_data, axis=1, copy=False)
            else:
                df = pd.DataFrame()
            fbs = encode_matrix_fbs(df, col_idx=df.columns, row_idx=None)

        return fbs

    def get_last_mod_time(self):
        try:
            lastmod = self.get_data_locator().lastmodtime()
        except RuntimeError:
            lastmod = None
        return lastmod
