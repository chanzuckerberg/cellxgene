"""
Helpers for user annotations / label_file parameter
"""
from os.path import exists, splitext
from os import remove, rename
import pandas as pd


def read_labels(fname):
    if exists(fname):
        return pd.read_csv(fname, dtype='category')
    else:
        return pd.DataFrame()


def write_labels(fname, df):
    rotate_fname(fname)
    df.to_csv(fname, index=False)


def rotate_fname(fname):
    """
    save N backups of file.
        fname -> fname-0
        fname-0 -> fname->1
        ...
        fname-(N-1) -> fname-N
    """

    def rotate(src, dst):
        if exists(src):
            if exists(dst):
                remove(dst)
            rename(src, dst)

    rotation_size = 9  # rotation size
    name, ext = splitext(fname)

    # rotate existing files
    for i in range(rotation_size - 1, 0, -1):
        src = f"{name}-{i}{ext}"
        tgt = f"{name}-{i+1}{ext}"
        rotate(src, tgt)

    tgt = f"{name}-1{ext}"
    rotate(fname, tgt)
