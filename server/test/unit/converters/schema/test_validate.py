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
