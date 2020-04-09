import unittest
from server.common.app_config import AppConfig

# NOTE, there are more tests that should be written for AppConfig.
# this is just a start.


class AppConfigTest(unittest.TestCase):
    def test_update(self):
        c = AppConfig()
        c.update(server__verbose=True, multi_dataset__dataroot="datadir")

        v = c.changes_from_default()
        self.assertCountEqual(v, [("server__verbose", True, False), ("multi_dataset__dataroot", "datadir", None)]),
