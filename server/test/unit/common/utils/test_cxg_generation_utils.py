import unittest
from os import popen, path, mkdir
from shutil import rmtree
from uuid import uuid4

import tiledb

from server.common.utils.cxg_generation_utils import convert_dictionary_to_cxg_group

PROJECT_ROOT = popen("git rev-parse --show-toplevel").read().strip()


class TestCxgGenerationUtils(unittest.TestCase):
    def setUp(self):
        self.testing_cxg_temp_directory = f"{PROJECT_ROOT}/server/test/fixtures/{uuid4()}"
        mkdir(self.testing_cxg_temp_directory)

    def tearDown(self):
        if path.isdir(self.testing_cxg_temp_directory):
            rmtree(self.testing_cxg_temp_directory)

    def test__convert_dictionary_to_cxg_group__writes_successfully(self):
        random_dictionary = {"cookies": "chocolate_chip", "brownies": "chocolate", "cake": "double chocolate"}
        dictionary_name = "favorite_desserts"
        expected_array_directory = f"{self.testing_cxg_temp_directory}/{dictionary_name}"

        convert_dictionary_to_cxg_group(self.testing_cxg_temp_directory, random_dictionary,
                                        group_metadata_name=dictionary_name)

        array = tiledb.open(expected_array_directory)
        actual_stored_metadata = dict(array.meta.items())

        self.assertTrue(path.isdir(expected_array_directory))
        self.assertTrue(isinstance(array, tiledb.DenseArray))
        self.assertEqual(random_dictionary, actual_stored_metadata)
