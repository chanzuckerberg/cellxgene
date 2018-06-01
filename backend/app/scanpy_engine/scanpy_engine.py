from os.path import join

import scanpy.api as sc
import numpy as np

from ..util.schema_parse import parse_schema

class ScanpyEngine():

	def __init__(self, dataloc):
		self.ADATA = sc.read(join(dataloc, "data.h5ad"))
		self.cell_count = self._cell_count
		self.schema = parse_schema(join(dataloc, "data_schema.json"))

	def _cell_count(self):
		return len(self.ADATA.obs.index)

	def all_cells(self):
		return self.ADATA.obs.index.tolist()

	def all_genes(self):
		return self.ADATA.var.index.tolist()

	def gene_count(self):
		return len(self.ADATA.var.index)

	def filter_cells(self, filter):
		cell_idx = np.ones((self.cell_count(),), dtype=bool)
		for key, value in filter.items():
			if value["variable_type"] == "categorical":
				key_idx = np.in1d(getattr(self.ADATA.obs, key), value["query"])
				cell_idx = np.logical_and(cell_idx, key_idx)
			else:
				min_ = value["query"]["min"]
				max_ = value["query"]["max"]
				if min_:
					key_idx = np.array((getattr(self.ADATA.obs, key) >= min_).data)
					cell_idx = np.logical_and(cell_idx, key_idx)
				if max_:
					key_idx = np.array((getattr(self.ADATA.obs, key) <= min_).data)
					cell_idx = np.logical_and(cell_idx, key_idx)
		return self.ADATA[cell_idx, :]

	@staticmethod
	def metadata_ranges(data, schema):
		metadata_ranges = {}
		for field in schema:
			if schema[field]["variabletype"] == "categorical":
				group_by = field
				if group_by == "CellName":
					group_by = 'index'
				metadata_ranges[field] = {"options": data.obs.groupby(group_by).size().to_dict()}
			else:
				metadata_ranges[field] = {
					"range": {
						"min": data.obs[field].min(),
						"max": data.obs[field].max()
					}
				}
		return metadata_ranges

	@staticmethod
	def metadata(data):
		cell_ids = []
		metadata = data.obs.to_dict(orient="records")
		# Do i have to loop twice?
		for idx, cell_name in enumerate(data.obs.index):
			metadata[idx]["CellName"] = cell_name
			cell_ids.append(cell_name)
		return metadata, cell_ids

	@staticmethod
	# TODO cache this
	# TODO accept n-dim versions too
	# TODO allow optional kw params to function
	def create_graph(data, graph_method="umap"):
		# Run the graph method
		getattr(sc.tl, graph_method)(data)
		graph = data.obsm["X_{graph_method}".format(graph_method=graph_method)]
		normalized_graph = (graph - graph.min()) / (graph.max() - graph.min())
		return np.hstack((data.obs.index.values.reshape(len(data.obs.index), 1), normalized_graph)).tolist()





