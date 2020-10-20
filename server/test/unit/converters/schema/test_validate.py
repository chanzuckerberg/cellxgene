import json
import unittest

import pandas as pd

from server.converters.schema import validate


class TestFieldValidation(unittest.TestCase):

    def test_validate_stringified_list_of_dicts(self):

        good = json.dumps([{"a": 1}, {2: "x", "z": "y"}])
        not_stringified = [{"a": 1}, {2: "x", "z": "y"}]
        not_a_list = json.dumps({"bad": "dict"})
        not_json = "oh hey!"

        self.assertTrue(validate._validate_stringified_list_of_dicts(good))

        self.assertFalse(validate._validate_stringified_list_of_dicts(not_stringified))
        self.assertFalse(validate._validate_stringified_list_of_dicts(not_a_list))
        self.assertFalse(validate._validate_stringified_list_of_dicts(not_json))

    def test_validate_human_readable_string(self):

        good = "oh hey!"
        curie = "EFO:0001"
        ensg = "ENSG000001234"
        enst = "ENST000005678"

        self.assertTrue(validate._validate_human_readable_string(good))

        self.assertFalse(validate._validate_human_readable_string(curie))
        self.assertFalse(validate._validate_human_readable_string(ensg))
        self.assertFalse(validate._validate_human_readable_string(enst))

    def test_validate_curie(self):

        self.assertTrue(validate._validate_curie("UBERON:00001", ["UBERON", "EFO"]))
        self.assertTrue(validate._validate_curie("HsapDv:00002", ["HsapDv"]))

        self.assertFalse(validate._validate_curie("HsapDv:00002", ["UBERON", "EFO"]))
        self.assertFalse(validate._validate_curie("EFO:00002 (organoid)", ["UBERON", "EFO"]))
        self.assertFalse(validate._validate_curie("EFO:00002 extra", ["UBERON", "EFO"]))
        self.assertFalse(validate._validate_curie("UBERON:ABCD", ["UBERON", "EFO"]))
        self.assertFalse(validate._validate_curie("Uberon:00002", ["UBERON", "EFO"]))
        self.assertFalse(validate._validate_curie("UBERON:", ["UBERON", "EFO"]))
        self.assertFalse(validate._validate_curie("UBERON", ["UBERON", "EFO"]))

    def test_validate_suffixed_curie(self):

        self.assertTrue(validate._validate_suffixed_curie("EFO:00001", ["UBERON", "EFO"]))
        self.assertTrue(validate._validate_suffixed_curie("UBERON:00001 (cell culture)", ["UBERON", "EFO"]))

        self.assertFalse(validate._validate_suffixed_curie("HsapDv:00002 (organoid)", ["UBERON", "EFO"]))
        self.assertFalse(validate._validate_suffixed_curie("HsapDv:00002(organoid)", ["UBERON", "EFO"]))
        self.assertFalse(validate._validate_suffixed_curie("HsapDv:00002", ["UBERON", "EFO"]))
        self.assertFalse(validate._validate_suffixed_curie("EFO:00002 extra", ["UBERON", "EFO"]))
        self.assertFalse(validate._validate_suffixed_curie("UBERON:ABCD", ["UBERON", "EFO"]))
        self.assertFalse(validate._validate_suffixed_curie("Uberon:00002", ["UBERON", "EFO"]))
        self.assertFalse(validate._validate_suffixed_curie("UBERON:", ["UBERON", "EFO"]))
        self.assertFalse(validate._validate_suffixed_curie("UBERON", ["UBERON", "EFO"]))


class TestColumnValidation(unittest.TestCase):

    def test_validate_unique(self):
        unique = pd.DataFrame([["abc", "def"], ["ghi", "jkl"], ["mnop", "qrs"]],
                              index=["X", "Y", "Z"], columns=["col1", "col2"])
        duped = pd.DataFrame([["abc", "def"], ["ghi", "qrs"], ["abc", "qrs"]],
                             index=["X", "Y", "X"], columns=["col1", "col2"])

        schema_def = {"unique": True}

        errors = validate._validate_column(unique.index, "index", "unique_df", schema_def)
        self.assertFalse(errors)

        errors = validate._validate_column(duped.index, "index", "duped_df", schema_def)
        self.assertEqual(len(errors), 1)
        self.assertIn("is not unique", errors[0])

        errors = validate._validate_column(unique["col1"], "col1", "unique_df", schema_def)
        self.assertFalse(errors)

        errors = validate._validate_column(duped["col1"], "col1", "duped_df", schema_def)
        self.assertEqual(len(errors), 1)
        self.assertIn("is not unique", errors[0])

        schema_def = {"unique": False}
        errors = validate._validate_column(duped["col1"], "col1", "duped_df", schema_def)
        self.assertFalse(errors)

    def test_validate_nullable(self):
        non_null = pd.DataFrame([["abc", "def"], ["ghi", "jkl"], ["mnop", "qrs"]],
                                index=["X", "Y", "Z"], columns=["col1", "col2"])
        has_null = pd.DataFrame([["abc", "", None], ["ghi", "jkl", 1], ["mnop", "qrs", 2]],
                                index=["X", "Y", "Z"], columns=["col1", "col2", "col3"])

        schema_def = {"nullable": False}
        errors = validate._validate_column(non_null["col1"], "col1", "nonnull_df", schema_def)
        self.assertFalse(errors)
        errors = validate._validate_column(has_null["col1"], "col1", "hasnull_df", schema_def)
        self.assertFalse(errors)

        errors = validate._validate_column(has_null["col2"], "col2", "hasnull_df", schema_def)
        self.assertEqual(len(errors), 1)
        self.assertIn("contains empty values", errors[0])
        errors = validate._validate_column(has_null["col3"], "col3", "hasnull_df", schema_def)
        self.assertEqual(len(errors), 1)
        self.assertIn("contains empty values", errors[0])

        schema_def = {"nullable": True}
        errors = validate._validate_column(has_null["col2"], "col2", "hasnull_df", schema_def)
        self.assertFalse(errors)
