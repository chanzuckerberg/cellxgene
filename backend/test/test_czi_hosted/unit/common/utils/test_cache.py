import concurrent
import unittest
from time import sleep
from unittest.mock import patch, MagicMock, Mock

from backend.czi_hosted.data_common.cache import CacheItem, CacheManager


class CacheItemTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.test_data = {"some": "data"}
        cls.item_one = CacheItem()
        cls.item_one.data = cls.test_data
        cls.item_two = CacheItem()
        cls.item_three = CacheItem()

    def create_data(self, cache_key, **create_data_args):
        data = {cache_key: cache_key}
        for key, value in create_data_args.items():
            data[key] = value
        return data

    def create_data_with_pause(self, cache_key, create_date_args):
        sleep(3)
        return self.create_data(cache_key=cache_key, create_date_args=create_date_args)

    def test_get_item_returns_data(self):
        data = self.item_one.get("key")
        self.assertEqual(self.test_data, data)

    def test_get_item_calls_create_data_lambda_with_all_args(self):
        data = self.item_two.get("cache_key", self.create_data)
        self.assertEqual({"cache_key":"cache_key"}, data)

        data = self.item_three.get("cache_key", self.create_data, {"a":"b", "c": "d"})
        self.assertEqual({"cache_key":"cache_key", "a":"b", "c": "d"}, data)

    def test_get_item_returns_none_if_no_data_and_no_create_data_lambda(self):
        cache_item = CacheItem()
        data = cache_item.get("cache_key")
        self.assertIsNone(data)

    def test_get_item_holds_write_lock_when_retrieving_data(self):
        cache_item = CacheItem()
        self.assertEqual(cache_item.data_lock.num_r, 0)
        with concurrent.futures.ThreadPoolExecutor(max_workers=2) as executor:
            executor.submit(cache_item.get, "cache_key", self.create_data_with_pause)
        try:
            cache_item.data_lock.w_acquire()
        self.assertEqual(cache_item.data_lock.num_r, 1)
        sleep(6)
        # self.assertEqual(cache_item.data_lock.num_r, 0)
        # import pdb
        # pdb.set_trace()
        cache_item_two = CacheItem()
        mock_create_data = Mock(spec=self.create_data)
        mock_create_data_with_pause = Mock(spec=self.create_data_with_pause)
        with concurrent.futures.ThreadPoolExecutor(max_workers=2) as executor:
            executor.submit(cache_item_two.get, "cache_key_two", mock_create_data_with_pause)
            executor.submit(cache_item_two.get, "cache_key_two", mock_create_data)
        self.assertEqual(mock_create_data_with_pause.call_count, 1)
        # the data should come in while this call is waiting for the lock so it shouldnt be necessary to call create_data
        self.assertEqual(mock_create_data.call_count, 0)
        import pdb
        pdb.set_trace()


    def test_write_lock_demoted_to_read_lock_when_data_acquired(self):
        pass

    def test_write_lock_demoted_when_exception_thrown_retrieving_data(self):
        pass

    def test_data_not_deleted_when_lock_held(self):
        pass

    def test_deletion_idempotent(self):
        cache_item = CacheItem()
        cache_item.get("cache_key", self.create_data)
        cache_item.attempt_delete()
        # this should not throw an error
        cache_item.attempt_delete()

class CacheManagerTest(unittest.TestCase):
    @classmethod
    def setUp(cls) -> None:
        cls.cache_manager = CacheManager(max_cached=3, timelimit_s=2)

    def test_cache_manager_retrieves_cached_data(self):
        get_data = MagicMock()
        get_data.return_value = {"cache_key": "cache_key"}
        with self.cache_manager.data_adaptor("key_one", lambda x: {x: x}) as cache_item:
            self.assertIsNotNone(cache_item)
        with self.cache_manager.data_adaptor("key_one", get_data) as cache_item:
            self.assertIsNotNone(cache_item)
            self.assertEqual(get_data.call_count, 0)

    def test_cache_manager_deletes_old_datasets(self):
        get_data = MagicMock()
        get_data.return_value = {"cache_key": "cache_key"}
        with self.cache_manager.data_adaptor("key_one", lambda x: {x: x}) as cache_item:
            self.assertIsNotNone(cache_item)
        with self.cache_manager.data_adaptor("key_one", get_data) as cache_item:
            self.assertIsNotNone(cache_item)
            self.assertEqual(get_data.call_count, 0)
        sleep(3)
        with self.cache_manager.data_adaptor("key_one", get_data) as cache_item:
            self.assertIsNotNone(cache_item)
            self.assertEqual(get_data.call_count, 1)
        mock_attempt_delete = Mock(spec=self.cache_manager.data.get("key_one").cache_item.attempt_delete)
        self.cache_manager.data.get("key_one").cache_item.attempt_delete = mock_attempt_delete
        sleep(3)
        with self.cache_manager.data_adaptor("key_two", lambda x: {x: x}) as cache_item_two:
            self.assertIsNotNone(cache_item_two)
        self.assertEqual(mock_attempt_delete.call_count, 1)

    def test_cache_manager_handles_lru_deletion(self):
        with self.cache_manager.data_adaptor("key_one", lambda x: {x: x}) as cache_item:
            self.assertIsNotNone(cache_item)
        with self.cache_manager.data_adaptor("key_two", lambda x: {x: x}) as cache_item_two:
            self.assertIsNotNone(cache_item_two)
        with self.cache_manager.data_adaptor("key_three", lambda x: {x: x}) as cache_item_three:
            self.assertIsNotNone(cache_item_three)
        get_data = MagicMock()
        get_data.return_value = {"cache_key": "cache_key"}
        # check first item is still cached
        with self.cache_manager.data_adaptor("key_one", get_data) as cache_item:
            self.assertIsNotNone(cache_item)
            self.assertEqual(get_data.call_count, 0)
        mock_attempt_delete = Mock(spec=self.cache_manager.data.get("key_two").cache_item.attempt_delete)
        self.cache_manager.data.get("key_two").cache_item.attempt_delete = mock_attempt_delete
        # add an additional cache item, check that the lru item was deleted
        with self.cache_manager.data_adaptor("key_four", lambda x: {x: x}) as cache_item:
            self.assertIsNotNone(cache_item)
        self.assertEqual(mock_attempt_delete.call_count, 1)
        # check it isnt in the cache
        with self.cache_manager.data_adaptor("key_two", get_data) as cache_item:
            self.assertIsNotNone(cache_item)
            self.assertEqual(get_data.call_count, 1)

    def test_cache_manager_sets_and_releases_read_lock(self):
        with self.cache_manager.data_adaptor("key_one", lambda x: {x: x}) as cache_item:
            self.assertIsNotNone(cache_item)
            self.assertGreater(self.cache_manager.data.get("key_one").cache_item.data_lock.num_r, 0)
        self.assertEqual(self.cache_manager.data.get("key_one").cache_item.data_lock.num_r, 0)

    def test_cache_manager_updates_cache_item_info_on_access(self):
        with self.cache_manager.data_adaptor("key_one", lambda x: {x: x}) as cache_item:
            self.assertIsNotNone(cache_item)
        cache_item_last_access = self.cache_manager.data.get("key_one").last_access
        cache_item_access_count = self.cache_manager.data.get("key_one").num_access
        with self.cache_manager.data_adaptor("key_one", lambda x: {x: x}) as cache_item:
            self.assertIsNotNone(cache_item)
        updated_cache_item_info = self.cache_manager.data.get("key_one")
        self.assertGreater(updated_cache_item_info.last_access, cache_item_last_access)
        self.assertGreater(updated_cache_item_info.num_access, cache_item_access_count)


