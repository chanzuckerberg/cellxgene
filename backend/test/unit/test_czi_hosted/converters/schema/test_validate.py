import json
import unittest

import pandas as pd
import scanpy as sc

from backend.czi_hosted.converters.schema import validate

# PROJECT_ROOT = os.popen("git rev-parse --show-toplevel").read().strip()
from backend.czi_hosted.test import PROJECT_ROOT


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

    def test_human_readable(self):
        hr_df = pd.DataFrame(
            [["for you, a human", "UBERON:12345", "UBERON:1234 (thundercat)"],
             ["hope you're well", "bit of lungs", "brain"]],
            index=["ENSG00001", "ENSG00002"],
            columns=["good", "curie", "suffixed_curie"])

        schema_def = {"type": "human-readable string"}
        errors = validate._validate_column(hr_df["good"], "good", "hr", schema_def)
        self.assertFalse(errors)

        errors = validate._validate_column(hr_df["curie"], "curie", "hr", schema_def)
        self.assertEqual(len(errors), 1)
        self.assertIn("non-human-readable", errors[0])

        errors = validate._validate_column(hr_df["suffixed_curie"], "suffixed_curie", "hr", schema_def)
        self.assertEqual(len(errors), 1)
        self.assertIn("non-human-readable", errors[0])

        errors = validate._validate_column(hr_df.index, "ensg", "hr", schema_def)
        self.assertEqual(len(errors), 1)
        self.assertIn("non-human-readable", errors[0])

    def test_curie(self):

        curie_df = pd.DataFrame(
            [["EFO:00001", "HsapDv:00001 (cell culture)", "EFO:", "MONDO:0001 cell culture"],
             ["UBERON:00002", "HsapDv:00002 (organoid)", "EFO:12345", "MONDO:0002 (baba yaga)"],
             ["EFO:0000000005", "HsapDv:000004 (humanzee)", "EFO:000002", "MONDO:0004 (TMNT)"]],
            index=["X", "Y", "Z"],
            columns=["good", "good_suffix", "bad", "bad_suffix"])

        # Good
        schema_def = {"type": "curie", "prefixes": ["EFO", "UBERON"]}
        errors = validate._validate_column(curie_df["good"], "good", "curie_df", schema_def)
        self.assertFalse(errors)

        # Good suffix
        schema_def = {"type": "suffixed curie", "prefixes": ["HsapDv", "WHATEVER"]}
        errors = validate._validate_column(curie_df["good_suffix"], "good_suffix", "curie_df", schema_def)
        self.assertFalse(errors)

        # Bad prefix
        schema_def = {"type": "curie", "prefixes": ["EFO"]}
        errors = validate._validate_column(curie_df["good"], "good", "curie_df", schema_def)
        self.assertEqual(len(errors), 1)
        self.assertIn("invalid ontology", errors[0])
        self.assertIn("must be curies from one of these", errors[0])

        # Bad curies
        schema_def = {"type": "curie", "prefixes": ["EFO"]}
        errors = validate._validate_column(curie_df["bad"], "bad", "curie_df", schema_def)
        self.assertEqual(len(errors), 1)
        self.assertIn("invalid ontology", errors[0])

        # Bad suffixes
        schema_def = {"type": "suffixed curie", "prefixes": ["EFO"]}
        errors = validate._validate_column(curie_df["bad_suffix"], "bad_suffix", "curie_df", schema_def)
        self.assertEqual(len(errors), 1)
        self.assertIn("invalid ontology", errors[0])

    def test_enum(self):
        enum_df = pd.DataFrame(
            [["abc", "ghi"],
             ["def", "jkl"]],
            index=["X", "Y"],
            columns=["col1", "col2"])

        # All match
        schema_def = {"type": "string", "enum": ["abc", "def", "xyz"]}
        errors = validate._validate_column(enum_df["col1"], "col1", "enum_df", schema_def)
        self.assertFalse(errors)

        # Missing value
        schema_def = {"type": "string", "enum": ["abc", "xyz"]}
        errors = validate._validate_column(enum_df["col1"], "col1", "enum_df", schema_def)
        self.assertEqual(len(errors), 1)
        self.assertIn("unpermitted values", errors[0])


class TestDictValidations(unittest.TestCase):


    def test_key_presence(self):

        schema_def = {"keys": {"abc": None, "def": None}}

        dict_ = {"abc": "123", "def": "456"}
        errors = validate._validate_dict(dict_, "d", schema_def)
        self.assertFalse(errors)

        # Missing keys are bad
        dict_ = {"abc": "123"}
        errors = validate._validate_dict(dict_, "d", schema_def)
        self.assertEqual(len(errors), 1)
        self.assertIn("missing key", errors[0])

        # Extra keys are okay
        dict_ = {"abc": "123", "def": "456", "xyz": "789"}
        errors = validate._validate_dict(dict_, "d", schema_def)
        self.assertFalse(errors)

        # Better not be empty come on
        dict_ = {}
        errors = validate._validate_dict(dict_, "d", schema_def)
        self.assertEqual(len(errors), 2)

    def test_nullable(self):

        schema_def = {"keys": {"abc": {"type": "string", "nullable": False},
                               "def": {"type": "string", "nullable": True}}}

        dict_ = {"abc": "xyz", "def": ""}
        errors = validate._validate_dict(dict_, "d", schema_def)
        self.assertFalse(errors)

        dict_ = {"abc": "", "def": ""}
        errors = validate._validate_dict(dict_, "d", schema_def)
        self.assertEqual(len(errors), 1)
        self.assertIn("empty value", errors[0])

    def test_recurse(self):

        schema_def = {
            "keys": {
                "subdict": {
                    "type": "dict",
                    "keys": {
                        "subdict_key1": None,
                        "subdict_key2": None
                    }
                },
                "ontology": {
                    "type": "curie",
                    "prefixes": ["ONTOLOGY"]
                },
                "blob": {
                    "type": "stringified list of dicts"
                }
            }
        }

        dict_ = {
            "subdict": {"subdict_key1": "any", "subdict_key2": "any"},
            "ontology": "ONTOLOGY:123456",
            "blob": json.dumps([{"abc": 123}, {"def": 456}])
        }
        errors = validate._validate_dict(dict_, "d", schema_def)
        self.assertFalse(errors)

        dict_ = {
            "subdict": {"subdict_key1": "any"},
            "ontology": "ONTOLOGY:123456",
            "blob": json.dumps([{"abc": 123}, {"def": 456}])
        }
        errors = validate._validate_dict(dict_, "d", schema_def)
        self.assertEqual(len(errors), 1)
        self.assertIn("missing key", errors[0])

        dict_ = {
            "subdict": {"subdict_key1": "any", "subdict_key2": "any"},
            "ontology": "oh no not an ontology term",
            "blob": json.dumps([{"abc": 123}, {"def": 456}])
        }
        errors = validate._validate_dict(dict_, "d", schema_def)
        self.assertEqual(len(errors), 1)
        self.assertIn("invalid ontology", errors[0])

        dict_ = {
            "subdict": {"subdict_key1": "any", "subdict_key2": "any"},
            "ontology": "ONTOLOGY:123456",
            "blob": [{"abc": 123}, {"def": 456}]
        }
        errors = validate._validate_dict(dict_, "d", schema_def)
        self.assertEqual(len(errors), 1)
        self.assertIn("JSON-encoded list of dicts", errors[0])

        # Multiple errors
        dict_ = {
            "subdict": {"subdict_key1": "any"},
            "ontology": "oh no not an ontology term",
            "blob": json.dumps([{"abc": 123}, {"def": 456}])
        }
        errors = validate._validate_dict(dict_, "d", schema_def)
        self.assertEqual(len(errors), 2)


class TestDataframeValidation(unittest.TestCase):

    def test_column_presence(self):
        df = pd.DataFrame(
            [["abc", "EFO:123"],
             ["def", "UBERON:456"]],
            columns=["hr_string", "ontology"],
            index=["X", "Y"]
        )

        schema_def = {
            "columns": {
                "hr_string": {"type": "human-readable string"},
                "ontology": {"type": "curie", "prefixes": ["EFO", "UBERON"]}
            }
        }
        errors = validate._validate_dataframe(df, "df", schema_def)
        self.assertFalse(errors)

        schema_def = {
            "columns": {
                "hr_string": {"type": "human-readable string"},
                "another_hr_string": {"type": "human-readable string"},
                "ontology": {"type": "curie", "prefixes": ["EFO", "UBERON"]}
            }
        }
        errors = validate._validate_dataframe(df, "df", schema_def)
        self.assertEqual(len(errors), 1)
        self.assertIn("missing column", errors[0])

        # Extra is okay
        df = pd.DataFrame(
            [["abc", "EFO:123", "extra"],
             ["def", "UBERON:456", "extra"]],
            columns=["hr_string", "ontology", "extra"],
            index=["X", "Y"]
        )
        schema_def = {
            "columns": {
                "hr_string": {"type": "human-readable string"},
                "ontology": {"type": "curie", "prefixes": ["EFO", "UBERON"]}
            }
        }
        errors = validate._validate_dataframe(df, "df", schema_def)
        self.assertFalse(errors)


    def test_index(self):
        df = pd.DataFrame(
            [["abc", "123"],
             ["def", "456"]],
            columns=["col1", "col2"],
            index=["ENSG0001", "ENSG0002"]
        )

        schema_def = {"index": {"unique": True}}
        errors = validate._validate_dataframe(df, "df", schema_def)
        self.assertFalse(errors)

        schema_def = {"index": {"type": "human-readable string"}}
        errors = validate._validate_dataframe(df, "df", schema_def)
        self.assertEqual(len(errors), 1)
        self.assertIn("non-human-readable", errors[0])

        df = pd.DataFrame(
            [["abc", "123"],
             ["def", "456"]],
            columns=["col1", "col2"],
            index=["ENSG0001", "ENSG0001"]
        )
        schema_def = {"index": {"unique": True}}
        errors = validate._validate_dataframe(df, "df", schema_def)
        self.assertEqual(len(errors), 1)
        self.assertIn("is not unique", errors[0])

    def test_recurse(self):

        df = pd.DataFrame(
            [["abc", "HsapDv:0001"],
             ["EFO:123", "UBERON:456"]],
            columns=["hr_string", "ontology"],
            index=["X", "Y"]
        )
        schema_def = {
            "columns": {
                "hr_string": {"type": "human-readable string"},
                "ontology": {"type": "curie", "prefixes": ["EFO", "UBERON"]}
            }
        }
        errors = validate._validate_dataframe(df, "df", schema_def)
        self.assertEqual(len(errors), 2)
        self.assertEqual(len([e for e in errors if "non-human-readable" in e]), 1)
        self.assertEqual(len([e for e in errors if "invalid ontology" in e]), 1)


class TestGetSchema(unittest.TestCase):

    def test_get_schema(self):
        self.assertIsInstance(validate.get_schema_definition("1.0.0"), dict)

        with self.assertRaises(ValueError):
            validate.get_schema_definition("10.1.5")


class TestValidate(unittest.TestCase):

    def setUp(self):
        self.source_h5ad_path = f"{PROJECT_ROOT}/backend/czi_hosted/test/fixtures/pbmc3k-CSC-gz.h5ad"

    def test_shallow(self):

        adata = sc.read_h5ad(self.source_h5ad_path)
        self.assertFalse(validate.validate_adata(adata, True))

        adata.uns["version"] = {
            "corpora_schema_version": "1.0.0",
            "corpora_encoding_version": "0.1.0"
        }
        self.assertTrue(validate.validate_adata(adata, True))

    def test_deep(self):
        adata = sc.read_h5ad(self.source_h5ad_path)
        self.assertFalse(validate.validate_adata(adata, False))

        adata.uns["version"] = {
            "corpora_schema_version": "1.0.0",
            "corpora_encoding_version": "0.1.0"
        }
        self.assertFalse(validate.validate_adata(adata, False))
