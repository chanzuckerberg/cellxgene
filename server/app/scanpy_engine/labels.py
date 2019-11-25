"""
Helpers for user annotations
"""
import os
import os.path
from datetime import datetime
import pandas as pd


def read_labels(fname):
    if fname is not None and os.path.exists(fname) and os.path.getsize(fname) > 0:
        return pd.read_csv(fname, dtype='category', index_col=0, header=0, comment='#')
    else:
        return pd.DataFrame()


def write_labels(fname, df, header=None, backup_dir=None):
    if backup_dir is not None:
        backup(fname, backup_dir)
        # rotate_fname(fname, backup_dir)
    if not df.empty:
        with open(fname, 'w', newline="") as f:
            if header is not None:
                f.write(header)
            df.to_csv(f)
    else:
        open(fname, 'w').close()


def backup(fname, backup_dir, max_backups=9):
    """
    save N backups of file to backup_dir.
        1. fname -> backup_dir/fname-TIME
        2. delete excess files in backup_dir
    """

    # Make sure there is work to do
    if not os.path.exists(fname):
        return

    # Ensure backup_dir exists
    if not os.path.exists(backup_dir):
        os.mkdir(backup_dir)

    # Save current file to backup_dir
    fname_base = os.path.basename(fname)
    fname_base_root, fname_base_ext = os.path.splitext(fname_base)
    # don't use ISO standard time format, as it contains characters illegal on some filesytems.
    nowish = datetime.now().strftime('%Y-%m-%dT%H-%M-%S')
    backup_fname = os.path.join(backup_dir, f"{fname_base_root}-{nowish}{fname_base_ext}")
    if os.path.exists(backup_fname):
        os.remove(backup_fname)
    os.rename(fname, backup_fname)

    # prune the backup_dir to max number of backup files, keeping the most recent backups
    backups = list(filter(lambda s: s.startswith(fname_base_root), os.listdir(backup_dir)))
    excess_count = len(backups) - max_backups
    if excess_count > 0:
        backups.sort()
        for bu in backups[0:excess_count]:
            os.remove(os.path.join(backup_dir, bu))
