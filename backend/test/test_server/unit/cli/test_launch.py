import filecmp
import os
import shutil
import unittest

import yaml

from backend.server.default_config import default_config
from backend.test import FIXTURES_ROOT


class CLIPLaunchTests(unittest.TestCase):
    tmp_dir = os.path.join(FIXTURES_ROOT, "dump_configs")

    @classmethod
    def setUpClass(cls) -> None:
        os.mkdir(cls.tmp_dir)

    @classmethod
    def tearDownClass(cls) -> None:
        shutil.rmtree(cls.tmp_dir)

    def test_dump_default_config(self):
        os.system(f"cellxgene launch --dump-default-config > {self.tmp_dir}/test_config_dump.txt")
        with open(f"{self.tmp_dir}/expected_config_dump.txt", "w") as expected_config:
            expected_config.write(yaml.dump(default_config))
        filecmp.cmp(f"{self.tmp_dir}/expected_config_dump.txt", f"{self.tmp_dir}/test_config_dump.txt")
