from enum import Enum
import threading
import time
from server.data_common.rwlock import RWLock
from server.common.errors import DatasetAccessError
from server.common.data_locator import DataLocator
from contextlib import contextmanager
from http import HTTPStatus


class MatrixDataCacheItem(object):
    """This class provides access and caching for a dataset.  The first time a dataset is accessed, it is
    opened and cached.  Later accesses use the cached version.   It may also be deleted by the
    MatrixDataCacheManager to make room for another dataset.  While a dataset is actively being used
    (during the lifetime of a api request), a reader lock is locked.  During that time, the dataset cannot
    be removed."""

    def __init__(self, loader):
        self.loader = loader
        self.data_adaptor = None
        self.data_lock = RWLock()

    def acquire_existing(self):
        """If the data_adaptor exists, take a read lock and return it, else return None"""
        self.data_lock.r_acquire()
        if self.data_adaptor:
            return self.data_adaptor

        self.data_lock.r_release()
        return None

    def acquire_and_open(self, app_config, dataset_config=None):
        """returns the data_adaptor if cached.  opens the data_adaptor if not.
        In either case, the a reader lock is taken.  Must call release when
        the data_adaptor is no longer needed"""
        self.data_lock.r_acquire()
        if self.data_adaptor:
            return self.data_adaptor
        self.data_lock.r_release()

        self.data_lock.w_acquire()
        # the data may have been loaded while waiting on the lock
        if not self.data_adaptor:
            try:
                self.loader.pre_load_validation()
                self.data_adaptor = self.loader.open(app_config, dataset_config)
            except Exception as e:
                # necessary to hold the reader lock after an exception, since
                # the release will occur when the context exits.
                self.data_lock.w_demote()
                raise DatasetAccessError(str(e))

        # demote the write lock to a read lock.
        self.data_lock.w_demote()
        return self.data_adaptor

    def release(self):
        """Release the reader lock"""
        self.data_lock.r_release()

    def delete(self):
        """Clear resources used by this dataset"""
        with self.data_lock.w_locked():
            if self.data_adaptor:
                self.data_adaptor.cleanup()
                self.data_adaptor = None

    def attempt_delete(self):
        """Delete, but only if the write lock can be immediately locked.  Return True if the delete happened"""
        if self.data_lock.w_acquire_non_blocking():
            if self.data_adaptor:
                try:
                    self.data_adaptor.cleanup()
                    self.data_adaptor = None
                except Exception:
                    # catch all exceptions to ensure the lock is released
                    pass

                self.data_lock.w_release()
                return True
        else:
            return False


class MatrixDataCacheInfo(object):
    def __init__(self, cache_item, timestamp):
        # The MatrixDataCacheItem in the cache
        self.cache_item = cache_item
        # The last time the cache_item was accessed
        self.last_access = timestamp
        # The number of times the cache_item was accessed (used for testing)
        self.num_access = 1


class MatrixDataCacheManager(object):
    """A class to manage the cached datasets.   This is intended to be used as a context manager
    for handling api requests.  When the context is created, the data_adator is either loaded or
    retrieved from a cache.  In either case, the reader lock is taken during this time, and release
    when the context ends.  This class currently implements a simple least recently used cache,
    which can delete a dataset from the cache to make room for a new one.

    This is the intended usage pattern:

           m = MatrixDataCacheManager(max_cached=..., timelimmit_s = ...)
           with m.data_adaptor(location, app_config) as data_adaptor:
               # use the data_adaptor for some operation
    """

    # FIXME:   If the number of active datasets exceeds the max_cached, then each request could
    # lead to a dataset being deleted and a new only being opened: the cache will get thrashed.
    # In this case, we may need to send back a 503 (Server Unavailable), or some other error message.

    # NOTE:  If the actual dataset is changed.  E.g. a new set of datafiles replaces an existing set,
    # then the cache will not react to this, however once the cache time limit is reached, the dataset
    # will automatically be refreshed.

    def __init__(self, max_cached, timelimit_s=None):
        # key is tuple(url_dataroot, location), value is a MatrixDataCacheInfo
        self.datasets = {}

        # lock to protect the datasets
        self.lock = threading.Lock()

        #  The number of datasets to cache.  When max_cached is reached, the least recently used
        #  cache is replaced with the newly requested one.
        #  TODO:  This is very simple.  This can be improved by taking into account how much space is actually
        #         taken by each dataset, instead of arbitrarily picking a max datasets to cache.
        self.max_cached = max_cached

        # items are automatically removed from the cache once this time limit is reached
        self.timelimit_s = timelimit_s

    @contextmanager
    def data_adaptor(self, url_dataroot, location, app_config):
        # create a loader for to this location if it does not already exist

        delete_adaptor = None
        data_adaptor = None
        cache_item = None

        key = (url_dataroot, location)
        with self.lock:
            self.evict_old_datasets()
            info = self.datasets.get(key)
            if info is not None:
                info.last_access = time.time()
                info.num_access += 1
                self.datasets[key] = info
                data_adaptor = info.cache_item.acquire_existing()
                cache_item = info.cache_item

            if data_adaptor is None:
                while True:
                    if len(self.datasets) < self.max_cached:
                        break

                    items = list(self.datasets.items())
                    items = sorted(items, key=lambda x: x[1].last_access)
                    # close the least recently used loader
                    oldest = items[0]
                    oldest_cache = oldest[1].cache_item
                    oldest_key = oldest[0]
                    del self.datasets[oldest_key]
                    delete_adaptor = oldest_cache

                loader = MatrixDataLoader(location, app_config=app_config)
                cache_item = MatrixDataCacheItem(loader)
                item = MatrixDataCacheInfo(cache_item, time.time())
                self.datasets[key] = item

        try:
            assert cache_item
            if delete_adaptor:
                delete_adaptor.delete()
            if data_adaptor is None:
                dataset_config = app_config.get_dataset_config(url_dataroot)
                data_adaptor = cache_item.acquire_and_open(app_config, dataset_config)
            yield data_adaptor
        except DatasetAccessError:
            cache_item.release()
            with self.lock:
                del self.datasets[key]
                cache_item.delete()
            cache_item = None
            raise

        finally:
            if cache_item:
                cache_item.release()

    def evict_old_datasets(self):
        # must be called with the lock held
        if self.timelimit_s is None:
            return

        now = time.time()
        to_del = []
        for key, info in self.datasets.items():
            if (now - info.last_access) > self.timelimit_s:
                # remove the data_cache when if it has been in the cache too long
                to_del.append((key, info))

        for key, info in to_del:
            # try and get the write_lock for the dataset.
            # if this returns false, it means the dataset is being used, and should
            # not be removed.
            if info.cache_item.attempt_delete():
                del self.datasets[key]


class MatrixDataType(Enum):
    H5AD = "h5ad"
    CXG = "cxg"
    UNKNOWN = "unknown"


class MatrixDataLoader(object):
    def __init__(self, location, matrix_data_type=None, app_config=None):
        """ location can be a string or DataLocator """
        region_name = None if app_config is None else app_config.server_config.data_locator__s3__region_name
        self.location = DataLocator(location, region_name=region_name)
        if not self.location.exists():
            raise DatasetAccessError("Dataset does not exist.", HTTPStatus.NOT_FOUND)

        # matrix_data_type is an enum value of type MatrixDataType
        self.matrix_data_type = matrix_data_type
        # matrix_type is a DataAdaptor type, which corresponds to the matrix_data_type
        self.matrix_type = None

        if matrix_data_type is None:
            self.matrix_data_type = self.__matrix_data_type()

        if not self.__matrix_data_type_allowed(app_config):
            raise DatasetAccessError("Dataset does not have an allowed type.")

        if self.matrix_data_type == MatrixDataType.CXG:
            from server.data_cxg.cxg_adaptor import CxgAdaptor

            self.matrix_type = CxgAdaptor

    def __matrix_data_type(self):
        if self.location.path.endswith(".h5ad"):
            return MatrixDataType.H5AD
        elif ".cxg" in self.location.path:
            return MatrixDataType.CXG
        else:
            return MatrixDataType.UNKNOWN

    def __matrix_data_type_allowed(self, app_config):
        if self.matrix_data_type == MatrixDataType.UNKNOWN or self.matrix_data_type == MatrixDataType.H5AD:
            return False

        if not app_config:
            return True
        if not app_config.is_multi_dataset():
            return True
        if not len(app_config.server_config.multi_dataset__allowed_matrix_types):
            return True

        for val in app_config.server_config.multi_dataset__allowed_matrix_types:
            try:
                if self.matrix_data_type == MatrixDataType(val):
                    return True
            except ValueError:
                # Check case where multi_dataset_allowed_matrix_type does not have a
                # valid MatrixDataType value.  TODO:  Add a feature to check
                # the AppConfig for errors on startup
                return False

        return False

    def pre_load_validation(self):
        if self.matrix_data_type == MatrixDataType.UNKNOWN:
            raise DatasetAccessError("Dataset does not have a recognized type: .h5ad or .cxg")
        self.matrix_type.pre_load_validation(self.location)

    def file_size(self):
        return self.matrix_type.file_size(self.location)

    def open(self, app_config, dataset_config=None):
        # create and return a DataAdaptor object
        return self.matrix_type.open(self.location, app_config, dataset_config)
