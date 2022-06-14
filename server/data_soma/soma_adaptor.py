import warnings
import tiledbsc
import numpy as np

import server.common.compute.diffexp_generic as diffexp_generic
import server.common.compute.estimate_distribution as estimate_distribution
from server.common.errors import PrepareError, DatasetAccessError
from server.common.colors import convert_soma_category_colors_to_cxg_category_colors
from server.common.fbs.matrix import encode_matrix_fbs
from server.common.constants import Axis, MAX_LAYOUTS, XApproximateDistribution
from server.common.utils.type_conversion_utils import get_schema_type_hint_of_array
from server.data_common.data_adaptor import DataAdaptor


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
            if data_locator.isfile():
                raise DatasetAccessError("is not a folder")

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
        # return dict(tiledbsc=str(tiledbsc.__version__))
        return dict(tiledbsc="0.0.0") # TODO: version information in tiledbsc library is broken

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
            elif not self._is_valid_layout(self.data.obsm[layout_name].df().to_numpy()):
                warnings.warn(f"Ignoring layout due to malformed shape or data type: {layout}")
            else:
                valid_layouts.append(layout)

        if len(valid_layouts) == 0:
            raise PrepareError("No valid layout data.")

        # cap layouts to MAX_LAYOUTS
        return valid_layouts[0:MAX_LAYOUTS]

    def get_embedding_array(self, ename, dims=2):
        """return an numpy array for the given pre-computed embedding name."""
        full_embedding = self.data.obsm[f"X_{ename}"].df().to_numpy()
        return full_embedding[:, 0:dims]

    def get_X_array(self, obs_mask=None, var_mask=None):
        """return the X array, possibly filtered by obs_mask or var_mask.
        the return type is either ndarray or scipy.sparse.spmatrix."""
        obs_ids = self.data.obs.df().index.to_numpy()
        var_ids = self.data.var.df().index.to_numpy()

        df = None
        obs_ids = obs_ids[obs_mask] if obs_mask is not None else None
        var_ids = var_ids[var_mask] if var_mask is not None else None

        # df = self.data.X.data.df(obs_ids, var_ids).reset_index() -> TODO: this tiledbsc function does not seem to work
        df = self.data.X.data.df().reset_index()
        df = df.loc[df["obs_id"].isin(obs_ids)] if obs_ids is not None else df
        df = df.loc[df["var_id"].isin(var_ids)] if var_ids is not None else df
        df = df.set_index(["obs_id", "var_id"])["value"].unstack().reset_index()
        df = df.set_index("obs_id")
        return df.to_numpy()

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
        return self.data.X.data.shape()

    def query_var_array(self, term_var):
        return self.get_var_df()[term_var]

    def query_obs_array(self, term_var):
        return self.get_obs_df()[term_var]

    def get_colors(self):
        return convert_soma_category_colors_to_cxg_category_colors(self.data)

    def get_obs_index(self):
        name = self.server_config.single_dataset__obs_names
        if name is None:
            return self.original_obs_index
        else:
            return self.data.obs.df()[name]

    def get_obs_columns(self):
        return self.data.obs.df().columns

    def get_obs_keys(self):
        # return list of keys
        return self.data.obs.keys()

    def get_var_keys(self):
        # return list of keys
        return self.data.var.keys()

    def get_obs_df(self):
        # get obs dataframe in the usual format, indexed by integer indices, not obs_id
        df = self.data.obs.df()
        df = df.rename_axis(self.get_schema()["annotations"]["obs"]["index"]).reset_index()
        return df

    def get_var_df(self):
        # get var dataframe in the usual format, indexed by integer indices, not var_id
        df = self.data.var.df()
        df = df.rename_axis(self.get_schema()["annotations"]["var"]["index"]).reset_index()
        return df

    def cleanup(self):
        pass

    def _create_schema(self):
        self.schema = {
            "dataframe": {
                "nObs": self.cell_count,
                "nVar": self.gene_count,
                **get_schema_type_hint_of_array(self.data.X.data.df().to_numpy()),
            },
            "annotations": {
                "obs": {"index": self.parameters.get("obs_names"), "columns": []},
                "var": {"index": self.parameters.get("var_names"), "columns": []},
            },
            "layout": {"obs": []},
        }
        for ax in Axis:
            curr_axis = None
            if str(ax) == "obs":
                curr_axis = self.get_obs_df()
            elif str(ax) == "var":
                curr_axis = self.get_var_df()
            else:
                raise RuntimeError("Invalid axis: %s." % str(ax))
            for soma in curr_axis:
                soma_schema = {"name": soma, "writable": False}
                soma_schema.update(get_schema_type_hint_of_array(curr_axis[soma]))
                self.schema["annotations"][ax]["columns"].append(soma_schema)

        for layout in self.get_embedding_names():
            layout_schema = {"name": layout, "type": "float32", "dims": [f"{layout}_0", f"{layout}_1"]}
            self.schema["layout"]["obs"].append(layout_schema)

    def get_schema(self):
        """
        Return current schema
        """
        return self.schema

    def annotation_to_fbs_matrix(self, axis, fields=None, labels=None):
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
                df.insert(0, self.schema["annotations"]["obs"]["index"], self.data.obs.df().index)
        else:
            df = self.data.var.df()
            df.insert(0, self.schema["annotations"]["var"]["index"], self.data.var.df().index)

        if fields is not None and len(fields) > 0:
            df = df[fields]

        df.index = range(df.shape[0])
        return encode_matrix_fbs(df, col_idx=df.columns)

    def compute_diffexp_ttest(self, maskA, maskB, top_n, lfc_cutoff):
        if top_n is None:
            top_n = self.dataset_config.diffexp__top_n
        if lfc_cutoff is None:
            lfc_cutoff = self.dataset_config.diffexp__lfc_cutoff
        return diffexp_generic.diffexp_ttest(self, maskA, maskB, top_n, lfc_cutoff)

    def _is_valid_layout(self, arr):
        """return True if this layout data is a valid array for front-end presentation:
        * ndarray, dtype float/int/uint
        * with shape (n_obs, >= 2)
        * with all values finite or NaN (no +Inf or -Inf)
        """
        is_valid = type(arr) == np.ndarray and arr.dtype.kind in "fiu"
        is_valid = is_valid and arr.shape[0] == self.data.obs.shape()[0] and arr.shape[1] >= 2
        is_valid = is_valid and not np.any(np.isinf(arr)) and not np.all(np.isnan(arr))
        return is_valid

    def _load_data(self, data_locator):
        try:
            # there is no guarantee data_locator indicates a local file. If we get a non-local object,
            # make a copy in tmp, and delete it after we load into memory.
            with data_locator.local_handle() as lh:
                self.data = tiledbsc.SOMA(lh)

        except ValueError:
            raise DatasetAccessError("Data must be in the SOMA format.")
        except MemoryError:
            raise DatasetAccessError("Out of memory - file is too large for available memory.")
        except Exception:
            raise DatasetAccessError(
                "Folder not found or is inaccessible. Folder must be in SOMA format."
                "Please check your input and try again."
            )

    def is_data_unique(self):
        obs_names = self.get_obs_keys()
        var_names = self.get_var_keys()
        return len(set(obs_names)) == len(obs_names) and len(set(var_names)) == len(var_names)

    def _validate_and_initialize(self):
        if not self.is_data_unique:
            raise KeyError("All annotation column names must be unique.")

        shape = self.get_shape()
        self._alias_annotation_names()
        self.cell_count = shape[0]
        self.gene_count = shape[1]
        self._create_schema()

        if self.dataset_config.X_approximate_distribution == "auto":
            self.get_X_approximate_distribution()
        else:
            self.X_approximate_distribution = self.dataset_config.X_approximate_distribution

        # heuristic
        n_values = shape[0] * shape[1]
        if n_values > 5e8:
            self.parameters.update({"diffexp_may_be_slow": True})

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
        self.original_obs_index = self.get_obs_columns()

        for (ax_name, var_name) in ((Axis.OBS, "obs"), (Axis.VAR, "var")):
            config_name = f"single_dataset__{var_name}_names"
            parameter_name = f"{var_name}_names"
            name = getattr(self.server_config, config_name)
            df_axis = getattr(self.data, str(ax_name)).df()
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

    def _create_unique_column_name(self, df, col_name_prefix):
        """given the columns of a dataframe, and a name prefix, return a column name which
        does not exist in the dataframe, AND which is prefixed by `prefix`

        The approach is to append a numeric suffix, starting at zero and increasing by
        one, until an unused name is found (eg, prefix_0, prefix_1, ...).
        """
        suffix = 0
        while f"{col_name_prefix}{suffix}" in df:
            suffix += 1
        return f"{col_name_prefix}{suffix}"
