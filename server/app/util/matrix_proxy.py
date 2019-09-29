import abc
from itertools import zip_longest
from copy import copy
import numpy as np

"""
cellxgene deals with a variety of matrix data types, many of which do
not support a consistent API.  This framework allows proxies to be created
to pave over some of this.  Most significantly, AnnData.X does not guarantee
that much of its API (eg. .X.T) will work.
"""

INT_TYPES = (int, np.integer)


class _ArrayProxyBase(abc.ABC):
    """
    Private base class for array or matrix proxy.  This summarizes
    the interface used by the rest of cellxgene.
    """
    @property
    @abc.abstractmethod
    def dtype(self):
        raise NotImplementedError()

    @property
    @abc.abstractmethod
    def ndim(self):
        raise NotImplementedError()

    @property
    @abc.abstractmethod
    def shape(self):
        raise NotImplementedError()

    @property
    @abc.abstractmethod
    def T(self):
        raise NotImplementedError()

    @abc.abstractmethod
    def __iter__(self):
        raise NotImplementedError()

    @abc.abstractmethod
    def __getitem__(self, args):
        raise NotImplementedError()

    @abc.abstractmethod
    def toarray():
        raise NotImplementedError()


class MatrixProxy(_ArrayProxyBase):
    """
    Abstract class - all interfaces we need, plus a factory method
    to create a proxy based upon actual matrix type.

    This class primarily provides the factory method and related support.
    All other functionality is delegated to subclasses.
    """

    """
    Registry of types to proxy class, where values are:
    * None: unsupported
    * True: self-supported
    * string: proxy class
    Sub-classes automatically register.
    """
    base_proxy_registry = {
        'pandas.core.frame.DataFrame': True,
        'numpy.ndarray': True,
        'scipy.sparse.csc.csc_matrix': True,
        'scipy.sparse.csr.csr_matrix': True,
    }
    proxy_registry = None
    last_cache_token = None

    @staticmethod
    def _register_subclasses(subclasses, registry):
        for c in subclasses:
            names = c.__supports__()
            for name in names:
                registry[name] = c
            MatrixProxy._register_subclasses(c.__subclasses__(), registry)

    @classmethod
    def build_proxy_registry(cls):
        if cls.proxy_registry and abc.get_cache_token() == cls.last_cache_token:
            return

        cls.last_cache_token = abc.get_cache_token()
        registry = copy(cls.base_proxy_registry)
        MatrixProxy._register_subclasses(cls.__subclasses__(), registry)
        cls.proxy_registry = registry

    @classmethod
    def create(cls, matrix):
        """
        Factory - call with a matrix and it will create a proxy if needed.
        If the type already supports the necessary API, it is just returned
        directly.
        """
        cls.build_proxy_registry()
        t = type(matrix)
        fqtn = t.__module__ + '.' + t.__name__
        proxy_cls = cls.proxy_registry.get(fqtn, None)
        if proxy_cls is None:
            raise Exception(f"Matrix format `{fqtn}` is unsupported by proxy.")
        if proxy_cls is True:
            return matrix
        return proxy_cls(matrix)

    def __init__(self, m):
        self.m = m

    @classmethod
    @abc.abstractmethod
    def __supports__(cls):
        raise NotImplementedError()

    @classmethod
    def ismatrixproxy(cls, m):
        return isinstance(m, _ArrayProxyBase)


class MatrixProxyView(MatrixProxy):
    """
    2D matrix view to a 2D matrix
    """
    def __init__(self, arg1, shape=None, index=(),
                 transposed=False, copy=False):
        if not copy:
            m = arg1
            super().__init__(m)

            if shape is None:
                shape = m.shape
            assert(len(shape) == 2)

            index = tuple(
                map(lambda s_i:
                    slice(0, s_i[0], 1) if s_i[1] is None else s_i[1],
                    zip_longest(shape, index))
            )

            self._shape = shape
            self._index = index
            self.transposed = transposed

        else:  # copy mode
            super().__init__(arg1.m)
            self._shape = arg1._shape
            self._index = arg1._index
            self.transposed = arg1.transposed

    @classmethod
    def create_array(cls, *args, **kwargs):
        """ override if you use a different 1D array proxy """
        return ArrayProxyView(*args, **kwargs)

    def copy(self):
        """ override if you need additional behaviors """
        return self.__class__(self, copy=True)

    @classmethod
    def __supports__(cls):
        return ()

    def _swap(self, x):
        return (x[1], x[0]) if self.transposed else x

    @property
    def dtype(self):
        return self.m.dtype

    @property
    def ndim(self):
        return len(self._shape)

    @property
    def shape(self):
        return self._swap(self._shape)

    @property
    def T(self):
        m = self.copy()
        m.transposed = not m.transposed
        return m

    def __iter__(self):
        M = self._swap(self._shape)[0]
        s = self._swap((self._index))[0]
        start, stop, step = s.indices(M)
        for d in range(start, stop, step):
            yield self[d]

    def __getitem__(self, args):
        """
        decompose into the indexing patterns we use, throw for the rest.
        Subclassses implement specialized access.
        """
        row, col = self._swap(_unpack_index(args, self.shape))
        M, N = self._shape
        allM, allN = self.m.shape

        if isinstance(row, INT_TYPES):
            row += self._index[0].start
        elif isinstance(row, slice):
            row = _slice_slice(self._index[0], allM, row, M)
        if isinstance(col, INT_TYPES):
            col += self._index[1].start
        elif isinstance(col, slice):
            col = _slice_slice(self._index[1], allN, col, N)

        if isinstance(row, INT_TYPES):
            if isinstance(col, INT_TYPES):
                return self._getitem_intXint(row, col)
            elif isinstance(col, slice):
                return self._getitem_intXslice(row, col)
        elif isinstance(row, slice):
            if isinstance(col, INT_TYPES):
                return self._getitem_sliceXint(row, col)
            elif isinstance(col, slice):
                return self._getitem_sliceXslice(row, col)

        raise IndexError("unsupported column index types")

    """
    These getitem signatures are separate so that they may be
    overridden by subclasses as necessary.  We don't do much
    with them by default other than the obvious sub-slicing.

    NOTE: these follow the numpy rules for dimensionality reduction
    when an integer index is specified.
    """
    def _getitem_intXint(self, row, col):
        return self.m[row, col]

    def _getitem_intXslice(self, row, col):
        shape = (_slice_length(col, self.m.shape[1]), )
        return self.__class__.create_array(self.m, shape=shape, index=(row, col))

    def _getitem_sliceXint(self, row, col):
        shape = (_slice_length(row, self.m.shape[0]), )
        return self.__class__.create_array(self.m, shape=shape, index=(row, col))

    def _getitem_sliceXslice(self, row, col):
        shape = (_slice_length(row, self.m.shape[0]),
                 _slice_length(col, self.m.shape[1]))
        return self.__class__(self.m, shape=shape, index=(row, col), transposed=self.transposed)

    def toarray(self):
        arr = self.m[self._index]
        if self.transposed:
            arr = arr.transpose()
        return arr


class ArrayProxyView(_ArrayProxyBase):
    """
    1D array view to a 2D matrix
    """
    def __init__(self, arg1, shape=None, index=None, copy=False):
        super().__init__()
        if not copy:
            m = arg1

            # one index MUST be an integer and the other MUST be a slice
            assert(len(index) == 2)
            assert(all(isinstance(idx, INT_TYPES + (slice, )) for idx in index))
            assert(isinstance(index[0], INT_TYPES) != isinstance(index[1], INT_TYPES))

            if shape is None:
                if isinstance(index[0], INT_TYPES):
                    shape = (m.shape[0], )
                else:
                    shape = (m.shape[1], )
            assert(len(shape) == 1)

            self._shape = shape
            self.m = m
            self._index = index
            self._vdim = 1 if isinstance(index[0], INT_TYPES) else 0

        else:
            self.m = arg1.m
            self._shape = arg1._shape
            self._index = arg1._index
            self.fixed = arg1.fixed

    def copy(self):
        return self.__class__(self, copy=True)

    @property
    def dtype(self):
        return self.m.dtype

    @property
    def ndim(self):
        return len(self._shape)

    @property
    def shape(self):
        return self._shape

    @property
    def T(self):
        return self.copy()

    def __iter__(self):
        _vdim = self._vdim
        M = self._shape[0]
        for d in range(*self._index[_vdim].indices(M)):
            yield self[d]

    def __getitem__(self, args):
        _vdim = self._vdim
        index = _unpack_index(args, self.shape)[0]
        M = self._shape[0]
        allM = self.m.shape[_vdim]

        if isinstance(index, INT_TYPES):
            index += self._index[_vdim].start
        elif isinstance(index, slice):
            index = _slice_slice(self._index[_vdim], allM, index, M)

        if _vdim == 0:
            row, col = index, self._index[1]
        else:
            row, col = self._index[0], index

        if isinstance(row, INT_TYPES):
            if isinstance(col, INT_TYPES):
                return self._getitem_intXint(row, col)
            elif isinstance(col, slice):
                return self._getitem_intXslice(row, col)
        elif isinstance(row, slice):
            assert(isinstance(col, INT_TYPES))
            return self._getitem_sliceXint(row, col)

        raise IndexError("unsupported column index types")

    def _getitem_intXint(self, row, col):
        return self.m[row, col]

    def _getitem_intXslice(self, row, col):
        shape = (_slice_length(col, self.m.shape[1]), )
        return self.__class__(self.m, shape=shape, index=(row, col))

    def _getitem_sliceXint(self, row, col):
        shape = (_slice_length(row, self.m.shape[0]), )
        return self.__class__(self.m, shape=shape, index=(row, col))

    def toarray(self):
        return self.m[self._index]


def _unpack_index(index, shape):
    if not isinstance(index, tuple):
        index = (index, )
    if len(shape) < len(index):
        raise IndexError("invalid index dimensionality - must be 2")

    unpacked = ()
    for shp, idx in zip_longest(shape, index):
        idx = slice(None) if idx is None else idx
        idx = _slice_defaults(idx, shp) if isinstance(idx, slice) else idx
        unpacked += (idx, )

    return unpacked


def _slice_slice(outer, outer_len, inner, inner_len):
    """
    slice a slice - we take advantage of Python 3 range's support
    for indexing.
    """
    assert(outer_len >= inner_len)
    outer_rng = range(*outer.indices(outer_len))
    rng = outer_rng[inner]
    start, stop, step = rng.start, rng.stop, rng.step
    if step < 0 and stop < 0:
        stop = None
    return slice(start, stop, step)


def _range_length(start, stop, step):
    """ return length of range """
    assert(step != 0)
    assert(start is not None and stop is not None and step is not None)
    if step > 0 and start < stop:
        return 1 + (stop - 1 - start) // step
    elif step < 0 and start > stop:
        return 1 + (start - 1 - stop) // -step
    else:
        return 0


def _slice_length(s, length):
    """ return slice length """
    return _range_length(*s.indices(length))


def _slice_defaults(s, length):
    """ apply slice defaulting conventions """
    assert(length >= 0)

    step = 1 if s.step is None else s.step

    if s.start is not None:
        start = s.start
        if start < 0:
            start += length
    else:
        start = 0 if step > 0 else (length - 1)

    if s.stop is not None:
        stop = s.stop
        if stop < 0:
            stop += length
    else:
        stop = length if step > 0 else -length - 1

    return slice(start, stop, step)
