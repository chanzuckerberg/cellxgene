import concurrent
import unittest
from time import sleep
from unittest.mock import patch, MagicMock, Mock

from backend.czi_hosted.data_common.cache import CacheItem, CacheManager


class CacheItemTest(unittest.TestCase):

    def create_data(self, cache_key, **create_data_args):
        data = {cache_key: cache_key}
        for key, value in create_data_args.items():
            data[key] = value
        return data

    def test_get_item_returns_data(self):
        cache_item = CacheItem()
        some_data = {"some": "data"}
        self.assertEqual(
            cache_item.get("key", lambda key: some_data), some_data
        )

    def test_get_item_calls_create_data_lambda_with_all_args(self):
        cache_item = CacheItem()
        data = cache_item.get("cache_key", self.create_data)
        self.assertEqual({"cache_key":"cache_key"}, data)

        cache_item = CacheItem()
        data = cache_item.get("cache_key", self.create_data, {"a":"b", "c": "d"})
        self.assertEqual({"cache_key":"cache_key", "a":"b", "c": "d"}, data)

    def test_get_item_returns_none_if_no_data_and_no_create_data_lambda(self):
        cache_item = CacheItem()
        data = cache_item.get("cache_key")
        self.assertIsNone(data)

    def test_get_locks_item(self):
        def assert_data_write_lock_held(cache_item, key):
            """ if called, write lock MUST be held, and read lock must NOT be held """
            self.assertTrue(cache_item.data_lock.w_lock.locked())
            self.assertFalse(cache_item.data_lock.num_r_lock.locked())
            self.assertEqual(cache_item.data_lock.num_r, 0)
            return "value"

        cache_item = CacheItem()

        # must hold write lock while creating item
        cache_item.get("key", lambda key: assert_data_write_lock_held(cache_item, key))

        # must hold read lock when item is returned
        self.assertEqual(cache_item.data_lock.num_r, 1)

        # and multiple readers must be allowed
        cache_item.get("key")
        self.assertEqual(cache_item.data_lock.num_r, 2)

        # and releasing all locks
        cache_item.release()
        cache_item.release()
        self.assertEqual(cache_item.data_lock.num_r, 0)
        self.assertFalse(cache_item.data_lock.num_r_lock.locked())
        self.assertFalse(cache_item.data_lock.w_lock.locked())

    def test_attempt_delete_deletes_unlocked_data(self):
        cache_item = CacheItem()
        some_data = {"some": "data"}
        data = cache_item.get("key", lambda key: some_data)
        self.assertIsNotNone(data)
        cache_item.release()

        deleted = cache_item.attempt_delete()
        self.assertTrue(deleted)

        # check data deleted
        data = cache_item.get("key")
        self.assertIsNone(data)

    def test_delete_calls_cleanup_and_clears_data(self):
        with self.subTest("Test delete call cleanup method on data adaptor"):
            class TestData:
                def cleanup(self):
                    pass

            some_data = TestData()
            some_data.cleanup = MagicMock()

            cache_item = CacheItem()
            data = cache_item.get("key", lambda key: some_data)
            self.assertIsNotNone(data)
            self.assertIsNotNone(cache_item.data)

            cache_item.release()
            cache_item.delete()
            self.assertEqual(some_data.cleanup.call_count, 1)
            self.assertIsNone(cache_item.data)

        with self.subTest("Test delete works without cleanup method"):
            cache_item = CacheItem()
            some_data = {"some": "data"}
            data = cache_item.get("key", lambda key: some_data)
            self.assertIsNotNone(data)
            self.assertIsNotNone(cache_item.data)

            cache_item.release()
            cache_item.delete()
            self.assertIsNone(cache_item.data)

    def test_release_function_releases_one_reader_lock(self):
        cache_item = CacheItem()
        cache_item.get("key", lambda key: {"some": "data"})

        self.assertEqual(cache_item.data_lock.num_r, 1)
        cache_item.get("key", lambda key: {"some": "data"})
        self.assertEqual(cache_item.data_lock.num_r, 2)

        cache_item.release()
        self.assertEqual(cache_item.data_lock.num_r, 1)

    def test_write_lock_demoted_to_read_lock_when_data_acquired(self):
        cache_item = CacheItem()
        some_data = {"some": "data"}
        cache_item.get("key", lambda key: some_data)
        self.assertFalse(cache_item.data_lock.d_lock.locked())
        self.assertEqual(cache_item.data_lock.num_r, 1)
        self.assertTrue(cache_item.data_lock.w_lock.locked())

    def test_write_lock_demoted_to_read_when_exception_thrown_retrieving_data(self):
        def raise_exception(cache_key):
            raise AssertionError

        cache_item = CacheItem()
        with self.assertRaises(AssertionError):
            cache_item.get("key", raise_exception)
        self.assertFalse(cache_item.data_lock.d_lock.locked())
        self.assertEqual(cache_item.data_lock.num_r, 1)
        self.assertTrue(cache_item.data_lock.w_lock.locked())

    def test_data_not_deleted_when_lock_held(self):
        cache_item = CacheItem()
        some_data = {"some": "data"}
        cache_item.get("key", lambda key: some_data)
        # check the lock
        self.assertEqual(cache_item.data_lock.num_r, 1)
        deleted = cache_item.attempt_delete()
        self.assertFalse(deleted)
        # check data still available
        data = cache_item.get("key")
        self.assertEqual(data, some_data)

        # drop locks and try deleting
        cache_item.release()
        cache_item.release()

        deleted = cache_item.attempt_delete()
        self.assertTrue(deleted)

        # check data deleted
        data = cache_item.get("key")
        self.assertIsNone(data)

    def test_deletion_idempotent(self):
        cache_item = CacheItem()
        cache_item.get("cache_key", self.create_data)
        cache_item.attempt_delete()
        # this should not throw an error
        cache_item.attempt_delete()

class CacheManagerTest(unittest.TestCase):

    def test_cache_manager_retrieves_cached_data(self):
        cache_manager = CacheManager(max_cached=3, timelimit_s=2)
        get_data = MagicMock()
        get_data.return_value = {"cache_key": "cache_key"}
        with cache_manager.data_adaptor("key_one", lambda x: {x: x}) as cache_item:
            self.assertIsNotNone(cache_item)
        with cache_manager.data_adaptor("key_one", get_data) as cache_item:
            self.assertIsNotNone(cache_item)
            self.assertEqual(get_data.call_count, 0)

    def test_cache_manager_deletes_old_datasets(self):
        cache_manager = CacheManager(max_cached=3, timelimit_s=2)
        get_data = MagicMock()
        get_data.return_value = {"cache_key": "cache_key"}
        with cache_manager.data_adaptor("key_one", lambda x: {x: x}) as cache_item:
            self.assertIsNotNone(cache_item)
        with cache_manager.data_adaptor("key_one", get_data) as cache_item:
            self.assertIsNotNone(cache_item)
            self.assertEqual(get_data.call_count, 0)
        sleep(3)
        with cache_manager.data_adaptor("key_one", get_data) as cache_item:
            self.assertIsNotNone(cache_item)
            self.assertEqual(get_data.call_count, 1)
        mock_attempt_delete = Mock(spec=cache_manager.data.get("key_one").cache_item.attempt_delete)
        cache_manager.data.get("key_one").cache_item.attempt_delete = mock_attempt_delete
        sleep(3)
        with cache_manager.data_adaptor("key_two", lambda x: {x: x}) as cache_item_two:
            self.assertIsNotNone(cache_item_two)
        self.assertEqual(mock_attempt_delete.call_count, 1)

    def test_cache_manager_handles_lru_deletion(self):
        cache_manager = CacheManager(max_cached=3, timelimit_s=2)

        with cache_manager.data_adaptor("key_one", lambda x: {x: x}) as cache_item:
            self.assertIsNotNone(cache_item)
        with cache_manager.data_adaptor("key_two", lambda x: {x: x}) as cache_item_two:
            self.assertIsNotNone(cache_item_two)
        with cache_manager.data_adaptor("key_three", lambda x: {x: x}) as cache_item_three:
            self.assertIsNotNone(cache_item_three)
        get_data = MagicMock()
        get_data.return_value = {"cache_key": "cache_key"}
        # check first item is still cached
        with cache_manager.data_adaptor("key_one", get_data) as cache_item:
            self.assertIsNotNone(cache_item)
            self.assertEqual(get_data.call_count, 0)
        mock_attempt_delete = Mock(spec=cache_manager.data.get("key_two").cache_item.attempt_delete)
        cache_manager.data.get("key_two").cache_item.attempt_delete = mock_attempt_delete
        # add an additional cache item, check that the lru item was deleted
        with cache_manager.data_adaptor("key_four", lambda x: {x: x}) as cache_item:
            self.assertIsNotNone(cache_item)
        self.assertEqual(mock_attempt_delete.call_count, 1)
        # check it isnt in the cache
        with cache_manager.data_adaptor("key_two", get_data) as cache_item:
            self.assertIsNotNone(cache_item)
            self.assertEqual(get_data.call_count, 1)

    def test_cache_manager_sets_and_releases_read_lock(self):
        cache_manager = CacheManager(max_cached=3, timelimit_s=2)

        with cache_manager.data_adaptor("key_one", lambda x: {x: x}) as cache_item:
            self.assertIsNotNone(cache_item)
            self.assertGreater(cache_manager.data.get("key_one").cache_item.data_lock.num_r, 0)
        self.assertEqual(cache_manager.data.get("key_one").cache_item.data_lock.num_r, 0)

    def test_cache_manager_updates_cache_item_info_on_access(self):
        cache_manager = CacheManager(max_cached=3, timelimit_s=2)

        with cache_manager.data_adaptor("key_one", lambda x: {x: x}) as cache_item:
            self.assertIsNotNone(cache_item)
        cache_item_last_access = cache_manager.data.get("key_one").last_access
        cache_item_access_count = cache_manager.data.get("key_one").num_access
        with cache_manager.data_adaptor("key_one", lambda x: {x: x}) as cache_item:
            self.assertIsNotNone(cache_item)
        updated_cache_item_info = cache_manager.data.get("key_one")
        self.assertGreater(updated_cache_item_info.last_access, cache_item_last_access)
        self.assertGreater(updated_cache_item_info.num_access, cache_item_access_count)


