import os
import unittest

import pandas as pd

from server.test import FIXTURES_ROOT
from server.converters.schema import gene_symbol


class TestHGNCSymbolChecker(unittest.TestCase):

    def setUp(self):
        self.test_hgnc_path = os.path.join(FIXTURES_ROOT, "hgnc_example.txt.gz")
        self.hgnc_checker = gene_symbol.HGNCSymbolChecker.from_hgnc_records(self.test_hgnc_path)

    def test_symbol_upgrade(self):
        self.assertEqual(self.hgnc_checker.upgrade_symbol("SEPT1"), "SEPTIN1")
        self.assertEqual(self.hgnc_checker.upgrade_symbol("ADRB2R"), "ADRB2")
        self.assertEqual(self.hgnc_checker.upgrade_symbol("BAR"), "ADRB2")
        self.assertEqual(self.hgnc_checker.upgrade_symbol("sept1"), "SEPTIN1")
        self.assertEqual(self.hgnc_checker.upgrade_symbol("AdRb2R"), "ADRB2")
        self.assertEqual(self.hgnc_checker.upgrade_symbol("bar"), "ADRB2")

        # DIFF6 is ambiguous so don't upgrade it
        self.assertEqual(self.hgnc_checker.upgrade_symbol("DIFF6"), "DIFF6")
        self.assertEqual(self.hgnc_checker.upgrade_symbol("diff6"), "diff6")

        # ARG1 is approved
        self.assertEqual(self.hgnc_checker.upgrade_symbol("ARG1"), "ARG1")
        self.assertEqual(self.hgnc_checker.upgrade_symbol("arg1"), "ARG1")

        # HAP1 is both approved and withdrawn
        self.assertEqual(self.hgnc_checker.upgrade_symbol("HAP1"), "HAP1")
        self.assertEqual(self.hgnc_checker.upgrade_symbol("hap1"), "HAP1")

        # Leave unknown symbols alone
        self.assertEqual(self.hgnc_checker.upgrade_symbol("NOTASYMBOL"), "NOTASYMBOL")
        self.assertEqual(self.hgnc_checker.upgrade_symbol("notasymbol"), "notasymbol")

    def test_check_symbol(self):
        self.assertEqual(self.hgnc_checker.check_symbol("SEPT1"), gene_symbol.SymbolStatus.UPGRADABLE)
        self.assertEqual(self.hgnc_checker.check_symbol("DIFF6"), gene_symbol.SymbolStatus.AMBIGUOUS)
        self.assertEqual(self.hgnc_checker.check_symbol("NOTASYMBOL"), gene_symbol.SymbolStatus.UNKNOWN)

        # HAP1 is one of the approved and withdrawn symbols
        self.assertEqual(self.hgnc_checker.check_symbol("HAP1"), gene_symbol.SymbolStatus.APPROVED)

    def test_upgrade_index(self):
        index = pd.Index(["SEPT1", "DIFF6", "NOTASYMBOL", "bar", "SEPTIN1"])
        var_df = pd.DataFrame([[0] * len(index)], index=index)
        upgraded_index = gene_symbol.get_upgraded_var_index(var_df, hgnc_path=self.test_hgnc_path)
        self.assertEqual(upgraded_index.tolist(), ["SEPTIN1", "DIFF6", "NOTASYMBOL", "ADRB2", "SEPTIN1"])
