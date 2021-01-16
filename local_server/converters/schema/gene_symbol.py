"""Helpers for converting and checking HGNC gene symbols."""

import argparse
import enum
import logging
import os
import re
import numpy as np
import pandas as pd


def get_upgraded_var_index(var, hgnc_path=None):
    """Given an anndata var dataframe, return a new index for the dataframe
    where human gene symbols have been upgraded to the current HGNC set.
    """

    if not hgnc_path:
        hgnc_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "hgnc_complete_set.txt.gz")

    hgnc_symbol_checker = HGNCSymbolChecker.from_hgnc_records(hgnc_path)

    return pd.Index([hgnc_symbol_checker.upgrade_symbol(s) for s in var.index])


class SymbolStatus(enum.Enum):
    """The status of a symbol in the HGNC database.

    APPROVED: Currently a valid symbol
    WITHDRAWN: A previously approved HGNC symbol for a gene that has since been shown
               not to exist _unless_ that symbol is also approved
    AMBIGUOUS: A symbol that is not approved but is an alias or previous symbol for
               multiple approved symbols
    UPGRADABLE: A symbol that is not approved but unambiguously maps to an approved
                symbol
    UNKNOWN: A symbol that does not appear in HGNC
    """

    APPROVED = 1
    WITHDRAWN = 2
    AMBIGUOUS = 3
    UPGRADABLE = 4
    UNKNOWN = 5


class HGNCSymbolChecker:
    """Handle checking and correcting HGNC symbols."""

    def __init__(self, approved_symbols, withdrawn_symbols, ambiguous_symbols, symbol_map):
        self.approved_symbols = approved_symbols
        self.withdrawn_symbols = withdrawn_symbols
        self.ambiguous_symbols = ambiguous_symbols
        self.symbol_map = symbol_map

    def print_symbol_map(self):
        """Print out a map from old symbol to new symbol."""

        for symbol_pair in self.symbol_map.items():
            print("\t".join(symbol_pair))

    def check_symbol(self, symbol):
        """See if a symbol if approved or something else."""
        if symbol in self.approved_symbols:
            return SymbolStatus.APPROVED

        if symbol in self.withdrawn_symbols:
            return SymbolStatus.WITHDRAWN

        if symbol in self.ambiguous_symbols:
            return SymbolStatus.AMBIGUOUS

        if symbol in self.symbol_map:
            return SymbolStatus.UPGRADABLE

        return SymbolStatus.UNKNOWN

    def upgrade_symbol(self, symbol):
        """Return the approved symbol for the given symbol.

        If the symbol cannot be upgraded, just return the original symbol.
        """

        fixed_symbol, stripped_symbol = format_symbol(symbol)

        if fixed_symbol in self.approved_symbols:
            return fixed_symbol
        elif fixed_symbol in self.symbol_map:
            return self.symbol_map[fixed_symbol]
        elif stripped_symbol in self.approved_symbols:
            return stripped_symbol
        elif stripped_symbol in self.symbol_map:
            return self.symbol_map[stripped_symbol]

        return symbol

    @classmethod
    def from_hgnc_records(cls, hgnc_dataset_path):
        """Parse a hgnc database download into a HGNCSymbolChecker object."""

        def all_symbols(record):
            """Get all the symbols associated with an HGNC record including previous, alias,
            and approved."""
            yield format_symbol(record["symbol"])[0]
            for symbol in alias_and_previous_symbols(record):
                yield symbol

        def alias_and_previous_symbols(record):
            """Get alias and previous symbols from an HGNC record."""
            for field in ("alias_symbol", "prev_symbol"):
                if record[field] is not np.nan:
                    for symbol in record[field].split("|"):
                        yield format_symbol(symbol)[0]
            # Sometimes something like HGNC:1234 appears in datasets, which we
            # want to fix as well.
            yield record["hgnc_id"]

        hgnc_records = pd.read_csv(hgnc_dataset_path, sep="\t", header=0, low_memory=False).to_dict("records")

        # Get all symbols that are currently approved.
        approved_symbols = set()
        for record in hgnc_records:
            if record["status"] == "Approved":
                approved_symbols.add(format_symbol(record["symbol"])[0])

        # Get all symbols that have been withdrawn
        withdrawn_symbols = set()
        for record in hgnc_records:
            if record["status"] == "Entry Withdrawn":
                for symbol in all_symbols(record):
                    withdrawn_symbols.add(symbol)

        # If a symbol is both approved and withdrawn, be optimistic and call it approved
        logging.warning(
            f"Some symbols are simulaneously withdrawn and approved\n"
            f"We will treat them at approved:\n"
            f"{withdrawn_symbols.intersection(approved_symbols)}"
        )
        withdrawn_symbols = withdrawn_symbols.difference(approved_symbols)

        # Now try to map from symbols that are not approved but are an alias or previous symbol for an approved symbol
        alias_previous_to_approved = {}
        ambiguous_symbols = set()

        for record in hgnc_records:
            if record["status"] == "Approved":

                # The approved symbol is what we'll map to
                approved_symbol = format_symbol(record["symbol"])[0]

                for symbol in alias_and_previous_symbols(record):

                    # If the alias or previous symbol is also an approved symbol,
                    # we'll just leave it alone
                    if symbol in approved_symbols:
                        continue

                    # If the alias or previous symbol maps to a different approved symbol, mark it as ambiguous
                    if symbol in alias_previous_to_approved and alias_previous_to_approved[symbol] != approved_symbol:
                        ambiguous_symbols.add(symbol)
                    else:
                        alias_previous_to_approved[symbol] = approved_symbol

        # Remove all the ambiguous symbols from the map
        for ambiguous_symbol in ambiguous_symbols:
            alias_previous_to_approved.pop(ambiguous_symbol)

        return HGNCSymbolChecker(approved_symbols, withdrawn_symbols, ambiguous_symbols, alias_previous_to_approved)


def format_symbol(symbol):
    """HGNC rules say symbols should all be upper case except for C#orf#. However, case is
    variable in both alias and previous symbols as well as in the symbols we get in
    submissions. So, upper case everything except for the one situation where mixed-case
    is allowed, which are the genes like C2orf157.

    Also, seurat and scanpy append ".1" or "-1" to duplicated gene names, and these altered
    names persist throughout the life of the object. They won't match against the HGNC database
    and we want to merge them, so we need to strip off the suffix and try matching again.

    This function takes a symbol and returns the symbol with the fixed case and also with the
    seurat/scanpy suffix stripped off.
    """

    match = re.match(r"^(C)(\d+)(orf)(\d+)$", symbol, re.IGNORECASE)

    if match:
        fixed_case = f"C{match.group(2)}orf{match.group(4)}"
    else:
        fixed_case = symbol.upper()

    suffix_stripped = re.sub(r"[\.\-]\d+$", "", fixed_case)

    return fixed_case, suffix_stripped


def main():
    """When called as main, parse a given hgnc download and print out a map from old to new
    symbol.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "hgnc_dataset", help="HGNC dataset tsv, available from www.genenames.org/download/statistics-and-files/"
    )
    args = parser.parse_args()

    hgnc_symbol_checker = HGNCSymbolChecker.from_hgnc_records(args.hgnc_dataset)

    hgnc_symbol_checker.print_symbol_map()


if __name__ == "__main__":
    main()
