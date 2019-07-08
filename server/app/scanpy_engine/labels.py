"""
Helpers for user annotations / label_file parameter
"""
from os.path import exists
import pandas as pd


def read_labels(fname):
    if exists(fname):
        return pd.read_csv(fname, dtype='category')
    else:
        return pd.DataFrame()


def write_labels(fname, df):
    """
    TODO / reminder:
    - need to rotate file if it already exists
    """
    df.to_csv(fname, index=False)
