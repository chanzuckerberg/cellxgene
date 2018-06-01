from abc import ABC, abstractmethod

class CXGDriver(metaclass=abc.ABCMeta):

	def __init__(self, data, schema=None, graph_method=None, diffexp_method=None):
		self.data = _load_data(data)

	@abstractmethod
	@staticmethod
	def _load_data(data):
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
		:return: iterator through cells
		"""
		pass

	@abstractmethod
	def metadata_ranges(self, cells_iterator):
		"""
		"""
		pass

	@abstractmethod
	def metadata(self, cells_iterator):
		"""
		"""
		pass

	@abstractmethod
	def create_graph(self, cells_iterator):
		"""
		"""
		pass


	@abstractmethod
	def diffexp(self, cells_iterator):
		"""
		"""
		pass
