from abc import ABCMeta, abstractmethod

class CXGDriver(metaclass=ABCMeta):

	def __init__(self, data, schema=None, graph_method=None, diffexp_method=None):
		self.data = self._load_data(data)

	@abstractmethod
	@staticmethod
	def _load_data(data):
		pass

	@abstractmethod
	def _load_or_infer_schema(data):
		pass

	@abstractmethod
	def _set_cell_ids(self):
		pass

	@abstractmethod
	def cells(self):
		pass

	@abstractmethod
	def cellids(self):
		pass

	@abstractmethod
	def genes(self):
		pass

	@abstractmethod
	def filter_cells(self, filter):
		"""
		Filter cells from data and return a subset of the data
		:param filter:
		:return: iterator through cell ids
		"""
		pass

	# Should this return the order of metadata fields as the first value?
	@abstractmethod
	def metadata(self, cells_iterator, fields=None):
		"""
		Generator for metadata. Gets the metadata values cell by cell and returns all value
		or only certain values if names is not None

		:param cells_iterator: from filter cells, iterator for cellids
		:param fields: list of keys for metadata to return, returns all metadata values if not set.
		:return: Iterator for cellid + list of cells metadata values  ex. [cell-id, val1, val2, val3]
		"""

		pass


	@abstractmethod
	def create_graph(self, cells_iterator):
		"""
		Computes a n-d layout for cells through dimensionality reduction.
		:param cells_iterator: from filter cells, iterator for cellids
		:return: Iterator for [cellid-1, pos1, pos2], [cellid-2, pos1, pos2]
		"""
		pass


	@abstractmethod
	def diffexp(self, cells_iterator_1, cells_iterator_2):
		"""
		Computes the top differentially expressed genes between two clusters

		:param cells_iterator_1: First set of cell ids
		:param cells_iterator_2: Second set of cell ids
		:return: Up in the air: I recommend [gene name, mean_expression_cells1, mean_expression_cells2, average_difference, statistic_value]
		"""
		pass
