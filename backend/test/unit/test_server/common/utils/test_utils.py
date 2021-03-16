import os
import shutil
import unittest

from backend.common.utils.utils import import_plugins
from backend.test.unit.test_server import PROJECT_ROOT, random_string


class TestPlugins(unittest.TestCase):
    """ Test plugin import functionality """

    plugins_dir = f"{PROJECT_ROOT}/backend/test/plugins"
    test_plugin_path = f"{plugins_dir}/foo.py"
    secret = random_string(8)

    @classmethod
    def setUpClass(cls) -> None:
        if not os.path.isdir(cls.plugins_dir):
            os.mkdir(cls.plugins_dir)
        with open(cls.test_plugin_path, "w") as fh:
            fh.write(f'SECRET = "{cls.secret}"\n')

    @classmethod
    def tearDownClass(cls) -> None:
        if os.path.isdir(cls.plugins_dir):
            shutil.rmtree(cls.plugins_dir)

    def test_import_plugins(self):
        self.assertTrue(os.path.isfile(self.test_plugin_path))
        loaded_modules = import_plugins("backend.test.plugins")
        # test that import plugins found the file
        self.assertEqual(["backend.test.plugins.foo"], [ele.__name__ for ele in loaded_modules])
        # test that the module was properly executed
        self.assertEqual(self.secret, loaded_modules[0].SECRET)
