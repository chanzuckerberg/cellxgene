import scanpy.api as sc
import numpy as np
import os

from ..util.schema_parse import parse_schema
from ..driver.driver import CXGDriver


class ScanpyEngine(CXGDriver):

	def __init__(self, data, schema=None, graph_method="umap", diffexp_method="ttest"):
		self.data = self._load_data(data)
		self.schema = self._load_or_infer_schema(data, schema)
		self._set_cell_ids()
		self.cell_count = self.data.shape[0]
		# TODO Do I need this?
		self.gene_count = self.data.shape[1]
		self.graph_method = graph_method
		self.diffexp_method = diffexp_method


	@staticmethod
	def _load_data(data):
		return sc.read(os.path.join(data, "data.h5ad"))

	@staticmethod
	def _load_or_infer_schema(data, schema):
		data_schema = None
		if not schema:
			pass
		else:
			data_schema = parse_schema(os.path.join(data,schema))
		return data_schema

	def _set_cell_ids(self):
		self.data.obs['cxg_cell_id'] = list(range(self.data.obs.shape[0]))
		self.data.obs["cell_name"] = list(self.data.obs.index)
		self.data.obs.set_index('cxg_cell_id', inplace=True)

	def cells(self):
		return list(self.data.obs.index)

	def cellids(self, cells_iterator=None):
		if cells_iterator:
			data = self.data.obs.iloc[[i for i in cells_iterator], :]
		else:
			data = self.data.obs
		return list(data.index)

	def genes(self):
		return self.data.var.index.tolist()

	def filter_cells(self, filter):
		"""
		Filter cells from data and return a subset of the data
		:param filter:
		:return: iterator through cell ids
		"""
		cell_idx = np.ones((self.cell_count,), dtype=bool)
		# TODO does this need to be a generator too?
		for key, value in filter.items():
			if value["variable_type"] == "categorical":
				key_idx = np.in1d(getattr(self.data.obs, key), value["query"])
				cell_idx = np.logical_and(cell_idx, key_idx)
			else:
				min_ = value["query"]["min"]
				max_ = value["query"]["max"]
				if min_:
					key_idx = np.array((getattr(self.data.obs, key) >= min_).data)
					cell_idx = np.logical_and(cell_idx, key_idx)
				if max_:
					key_idx = np.array((getattr(self.data.obs, key) <= min_).data)
					cell_idx = np.logical_and(cell_idx, key_idx)
		# If this is slow, could vectorize with logical array and then loop through that
		for idx in range(self.cell_count):
			if cell_idx[idx]:
				yield self.data.obs.index[idx]


	def metadata_ranges(self, cells_iterator=None):
		metadata_ranges = {}
		if cells_iterator:
			data = self.data.obs.iloc[[i for i in cells_iterator], :]
		else:
			data = self.data.obs
		for field in self.schema:
			if self.schema[field]["variabletype"] == "categorical":
				group_by = field
				if group_by == "CellName":
					group_by = 'cell_name'
				metadata_ranges[field] = {"options": data.groupby(group_by).size().to_dict()}
			else:
				metadata_ranges[field] = {
					"range": {
						"min": data[field].min(),
						"max": data[field].max()
					}
				}
		return metadata_ranges

	# Should this return the order of metadata fields as the first value?
	def metadata(self, cells_iterator, fields=None):
		"""
		Generator for metadata. Gets the metadata values cell by cell and returns all value
		or only certain values if names is not None

		:param cells_iterator: from filter cells, iterator for cellids
		:param fields: list of keys for metadata to return, returns all metadata values if not set.
		:return: Iterator for cellid + list of cells metadata values  ex. [cell-id, val1, val2, val3]
		"""
		if not fields:
			fields = self.data.obs.columns.tolist()
		for cell_id in cells_iterator:
			yield [cell_id] + self.data.obs.loc[cell_id, fields].tolist()


	def create_graph(self, cells_iterator):
		"""
		Computes a n-d layout for cells through dimensionality reduction.
		:param cells_iterator: from filter cells, iterator for cellids
		:return: Iterator for [cellid-1, pos1, pos2], [cellid-2, pos1, pos2]
		"""
		cell_ids = list(cells_iterator)
		getattr(sc.tl, self.graph_method)(self.data[self.data.obs.index.isin(cell_ids)])
		graph = self.data.obsm["X_{graph_method}".format(graph_method=self.graph_method)]
		normalized_graph = (graph - graph.min()) / (graph.max() - graph.min())
		for idx, cell_id in enumerate(cell_ids):
			yield [cell_id] + normalized_graph[idx].tolist()


	def diffexp(self, cells_iterator_1, cells_iterator_2):
		"""
		Computes the top differentially expressed genes between two clusters

		:param cells_iterator_1: First set of cell ids
		:param cells_iterator_2: Second set of cell ids
		:return: Up in the air: I recommend [gene name, mean_expression_cells1, mean_expression_cells2, average_difference, statistic_value]
		"""
		pass






