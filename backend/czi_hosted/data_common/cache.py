import threading
import time
from collections import Callable
from typing import Optional

import typing

from backend.common.errors import DatasetAccessError, DatasetNotFoundException
from contextlib import contextmanager

from backend.czi_hosted.data_common.rwlock import RWLock


class CacheItem(object):
    """This class provides a base class for caching information.  The first time information is accessed,
    it is located (via the lambda function passed to .get) and cached.  Later accesses use the cached version.
    It may also be deleted by the DataCacheManager to make room.  While data is actively being used
    (during the lifetime of an api request), a reader lock is locked.
    During that time, the data cannot be removed or updated."""

    def __init__(self):
        # self.loader = loader
        self.data_lock = RWLock()
        self.data = None


    def get(self, cache_key: str, create_data_lambda: typing.Optional[typing.Callable[[str], object]] = None, create_data_args: object={}):
        self.data_lock.r_acquire()
        if self.data:
            return self.data
        self.data_lock.r_release()
        if create_data_lambda is None:
            return None

        self.data_lock.w_acquire()
        # the data may have been loaded while waiting on the lock
        if not self.data:
            try:
                # self.data = self.loader.
                self.data = create_data_lambda(cache_key, **create_data_args)
            except Exception:
                # necessary to hold the reader lock after an exception, since
                # the release will occur when the context exits.
                self.data_lock.w_demote()
                raise

        # demote the write lock to a read lock.
        self.data_lock.w_demote()
        return self.data

    def release(self):
        """Release the reader lock"""
        self.data_lock.r_release()

    def delete(self):
        """Clear resources used by this cache item, if self.data has a cleanup method, call it, if not
        delete the cached information"""
        with self.data_lock.w_locked():
            if self.data:
                try:
                    self.data.cleanup()
                except AttributeError:
                    pass
                self.data = None

    def attempt_delete(self):
        """Delete, but only if the write lock can be immediately locked.  Return True if the delete happened"""
        if self.data_lock.w_acquire_non_blocking():
            if self.data:
                try:
                    self.data.cleanup()
                    self.data = None
                except AttributeError:
                    self.data = None
                except Exception:
                    # catch all exceptions to ensure the lock is released
                    pass

                self.data_lock.w_release()
                return True
        else:
            return False


class CacheItemInfo(object):
    """
    This class stores a pointer to and metadata about the cached item. When the data is accessd the
    update function is called to track the latest access
    """
    def __init__(self, cache_item: CacheItem(), timestamp):
        self.cache_item = cache_item
        self.last_access = timestamp
        self.num_access = 1

    def update_latest_cache_access(self):
        self.last_access = time.time()
        self.num_access +=1

class CacheManager(object):
    """An base class to manage cached data. This is intended to be used as a context manager
    for handling api requests.  When the context is created, the data_adaptor is either loaded or
    retrieved from a cache.  In either case, the reader lock is taken during this time, and release
    when the context ends.  This class currently implements a simple least recently used cache,
    which can delete data from the cache to make room for a new one.

    This is the intended usage pattern:

           cache_manager = DataCacheManager(max_cached=..., timelimmit_s = ...)
           with cache_manager.data_adaptor(location, app_config) as data_adaptor:
               # use the data_adaptor for some operation
    """

    # FIXME:   If the number of active datasets exceeds the max_cached, then each request could
    # lead to a dataset being deleted and a new only being opened: the cache will get thrashed.
    # In this case, we may need to send back a 503 (Server Unavailable), or some other error message.

    # NOTE:  If the actual dataset is changed.  E.g. a new set of datafiles replaces an existing set, or the location
    # of the data is updated then the cache will not react to this, however once the cache time limit is reached,
    # the data will automatically be refreshed.

    def __init__(self, max_cached: int, timelimit_s: int=None):
        self.data = {}

        # lock to protect the datasets
        self.lock = threading.Lock()

        #  The number of items to cache.  When max_cached is reached, the least recently used
        #  cache is replaced with the newly requested one.
        #  TODO:  This is very simple.  When this class is used to cache the actual datasets
        #   this can be improved by taking into account how much space is actually taken by each dataset,
        #   instead of arbitrarily picking a max number to cache.
        self.max_cached = max_cached

        # items are automatically removed from the cache once this time limit (in seconds) is reached
        self.timelimit_s = timelimit_s

    @contextmanager
    def data_adaptor(self, cache_key: str, create_data_lambda: Optional[Callable]=None, create_data_args: object={}):
        """"
        If no data is stored under the given cache_key, pass the create_data_lamba to the CacheItem.get method
        """

        delete_adaptor = None
        desired_data = None
        cache_item = None

        with self.lock:
            self.evict_old_data()
            cache_item_info = self.data.get(cache_key)
            if cache_item_info is not None:
                cache_item_info.update_latest_cache_access()
                cache_item = cache_item_info.cache_item

            if cache_item is None:
                delete_adaptor = self.evict_extra_data()
                cache_item = CacheItem()
                desired_data = cache_item.get(
                    cache_key=cache_key,
                    create_data_lambda=create_data_lambda,
                    create_data_args=create_data_args
                )
                cache_item_info = CacheItemInfo(cache_item, time.time())
                self.data[cache_key] = cache_item_info

        try:
            assert cache_item
            if delete_adaptor:
                delete_adaptor.delete()
            # TODO @mdunitz can this be removed?
            if desired_data is None:
                desired_data = cache_item.get(cache_key)

            yield desired_data
        except (DatasetAccessError, DatasetNotFoundException):
            cache_item.release()
            with self.lock:
                del self.data[cache_key]
                cache_item.delete()
            cache_item = None
            raise

        finally:
            if cache_item:
                cache_item.release()


    def evict_extra_data(self):
        delete_adaptor = None
        while True:
            if len(self.data) < self.max_cached:
                break

            items = list(self.data.items())
            items = sorted(items, key=lambda x: x[1].last_access)
            # close the least recently used loader
            oldest = items[0]
            oldest_cache = oldest[1].cache_item
            oldest_key = oldest[0]
            del self.data[oldest_key]
            delete_adaptor = oldest_cache
        return delete_adaptor

    def evict_old_data(self):
        # must be called with the lock held
        if self.timelimit_s is None:
            return

        now = time.time()
        to_del = []
        for key, info in self.data.items():
            if (now - info.last_access) > self.timelimit_s:
                # remove the data_cache when if it has been in the cache too long
                to_del.append((key, info))

        for key, info in to_del:
            # try and get the write_lock for the dataset.
            # if this returns false, it means the dataset is being used, and should
            # not be removed.
            if info.cache_item.attempt_delete():
                del self.data[key]
