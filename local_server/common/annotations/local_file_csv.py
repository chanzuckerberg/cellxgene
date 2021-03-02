import base64
import os
import re
import threading
from datetime import datetime
from hashlib import blake2b
import csv

import pandas as pd
from flask import session, has_request_context, current_app

from local_server import __version__ as cellxgene_version
from local_server.common.annotations.annotations import Annotations
from local_server.common.errors import AnnotationsError, ObsoleteRequest


class AnnotationsLocalFile(Annotations):
    CXG_ANNO_COLLECTION = "cxg_anno_collection"

    def __init__(self, config, output_dir, label_output_file, genesets_output_file):
        super().__init__(config)
        self.output_dir = output_dir
        self.label_output_file = label_output_file
        self.genesets_output_file = genesets_output_file
        # lock used to protect label file write ops
        self.label_lock = threading.RLock()
        self.genesets_lock = threading.RLock()

        # cache the most recent annotations.
        self.last_fname = None
        self.last_labels = None

        # txn ID - used to de-dup geneset writes
        self.last_geneset_tid = 0

    def is_safe_collection_name(self, name):
        """
        return true if this is a safe collection name
        this is ultra conservative. If we want to allow full legal file name syntax,
        we could look at modules like `pathvalidate`
        """
        if name is None:
            return False
        return re.match(r"^[\w\-]+$", name) is not None

    def set_collection(self, name):
        session[self.CXG_ANNO_COLLECTION] = name
        session.permanent = True

    def get_collection(self):
        if session is None:
            return None
        return session.get(self.CXG_ANNO_COLLECTION)

    def read_labels(self, data_adaptor):
        self.check_user_annotations_enabled()  # raises

        if has_request_context():
            if not current_app.auth.is_user_authenticated():
                return pd.DataFrame()

        fname = self._get_celllabels_filename(data_adaptor)
        with self.label_lock:
            if fname is not None and os.path.exists(fname) and os.path.getsize(fname) > 0:
                # returned the cached labels if possible, otherwise read them from the file
                if fname == self.last_fname:
                    return self.last_labels
                else:
                    labels = pd.read_csv(
                        fname, dtype="category", index_col=0, header=0, comment="#", keep_default_na=False
                    )
                    # update the cache
                    self.last_fname = fname
                    self.last_labels = labels
                    return labels
            else:
                return pd.DataFrame()

    def write_labels(self, df, data_adaptor):
        self.check_user_annotations_enabled()  # raises

        # update our internal state and save it.  Multi-threading often enabled,
        # so treat this as a critical section.
        with self.label_lock:
            lastmod = data_adaptor.get_last_mod_time()
            lastmodstr = "'unknown'" if lastmod is None else lastmod.isoformat(timespec="seconds")
            header = (
                f"# Annotations generated on {datetime.now().isoformat(timespec='seconds')} "
                f"using cellxgene version {cellxgene_version}\n"
                f"# Input data file was {data_adaptor.get_location()}, "
                f"which was last modified on {lastmodstr}\n"
            )

            fname = self._get_celllabels_filename(data_adaptor)
            self._backup(fname)
            if not df.empty:
                with open(fname, "w", newline="") as f:
                    if header is not None:
                        f.write(header)
                    df.to_csv(f)
            else:
                open(fname, "w").close()

            # update the cache
            self.last_fname = fname
            self.last_labels = df

    def read_genesets(self, data_adaptor, context=None):
        if has_request_context():
            if not current_app.auth.is_user_authenticated():
                return ({}, self.last_geneset_tid)

        fname = self._get_genesets_filename(data_adaptor)
        genesets = {}
        tid = None
        with self.genesets_lock:
            tid = self.last_geneset_tid  # inside the critical section
            if fname is not None and os.path.exists(fname) and os.path.getsize(fname) > 0:
                with open(fname, newline="") as f:
                    genesets = read_geneset_tidycsv(f, context)

        return (genesets, tid)

    def write_genesets(self, genesets, tid, data_adaptor):
        self.check_genesets_save_enabled()  # raises

        if type(tid) != int or tid < 0:
            raise ValueError("tid must be a positive integer")

        with self.genesets_lock:
            # skip if the request is stale
            if tid is not None:
                if tid <= self.last_geneset_tid:
                    raise ObsoleteRequest("TID is stale.")
                self.last_geneset_tid = tid

            lastmod = data_adaptor.get_last_mod_time()
            lastmodstr = "'unknown'" if lastmod is None else lastmod.isoformat(timespec="seconds")
            header = (
                f"# Geneset generated on {datetime.now().isoformat(timespec='seconds')} "
                f"using cellxgene version {cellxgene_version}\n"
                f"# Input data file was {data_adaptor.get_location()}, "
                f"which was last modified on {lastmodstr}\n"
            )

            fname = self._get_genesets_filename(data_adaptor)
            self._backup(fname)
            with open(fname, "w", newline="") as f:
                f.write(header)
                f.write(self.genesets_to_csv(genesets))

    def _get_userdata_idhash(self, data_adaptor):
        """
        Return a short hash that weakly identifies the user and dataset.
        Used to create safe annotations output file names.
        """
        uid = current_app.auth.get_user_id() or ""
        id = (uid + data_adaptor.get_location()).encode()
        idhash = base64.b32encode(blake2b(id, digest_size=5).digest()).decode("utf-8")
        return idhash

    def _get_output_dir(self):
        if self.output_dir:
            return self.output_dir

        output_file = self.label_output_file or self.genesets_output_file
        if output_file:
            return os.path.dirname(os.path.abspath(output_file))

        return os.getcwd()

    def _get_celllabels_filename(self, data_adaptor):
        """ return the current annotation file name """
        if self.label_output_file:
            return self.label_output_file

        return self._get_filename(data_adaptor, "celllabels")

    def _get_genesets_filename(self, data_adaptor):
        """ return the current genesets file name """
        if self.genesets_output_file:
            return self.genesets_output_file

        return self._get_filename(data_adaptor, "genesets")

    def _get_filename(self, data_adaptor, anno_name):
        # we need to generate a file name, which we can only do if we have a UID and collection name
        if session is None:
            raise AnnotationsError("unable to determine file name for annotations")

        collection = self.get_collection()
        if collection is None:
            return None

        if data_adaptor is None:
            raise AnnotationsError("unable to determine file name for annotations")

        idhash = self._get_userdata_idhash(data_adaptor)
        return os.path.join(self._get_output_dir(), f"{collection}-{anno_name}-{idhash}.csv")

    def _backup(self, fname, max_backups=9):
        """
        save N backups of file to backup_dir.
            1. fname -> backup_dir/fname-TIME
            2. delete excess files in backup_dir
        """
        root, ext = os.path.splitext(fname)
        backup_dir = f"{root}-backups"

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
        nowish = datetime.now().strftime("%Y-%m-%dT%H-%M-%S")
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

    def update_parameters(self, parameters, data_adaptor):
        params = {}
        params["annotations"] = self.user_annotations_enabled()
        params["annotations_genesets_readonly"] = not self.genesets_save_enabled()
        params["user_annotation_collection_name_enabled"] = True

        if self.ontology_data:
            params["annotations_cell_ontology_enabled"] = True
            params["annotations_cell_ontology_terms"] = self.ontology_data
        else:
            params["annotations_cell_ontology_enabled"] = False

        if self.label_output_file is not None:
            # user has hard-wired the name of the annotation cell label data collection
            fname = os.path.basename(self.label_output_file)
            collection_fname = os.path.splitext(fname)[0]
            params["annotations-data-collection-is-read-only"] = True
            params["annotations-data-collection-name"] = collection_fname

        elif session is not None:
            collection = self.get_collection()
            params["annotations-data-collection-is-read-only"] = False
            params["annotations-data-collection-name"] = collection

        if current_app.auth.is_user_authenticated():
            params["annotations-user-data-idhash"] = self._get_userdata_idhash(data_adaptor)

        parameters.update(params)


def read_geneset_tidycsv(f, context=None):
    """
    Read & parse the Tidy CSV format, applying validation checks for mandatory
    values, and de-duping rules.

    Format is a four-column CSV, with a mandatory header row, and optional "#" prefixed
    comments.  Format:

        geneset_name, geneset_description, gene_symbol, gene_description

    geneset_name and gene_symbol must be non-null; others are optional.

    Returns: a dictionary of the shape (values in angle-brackets vary):

        {
            <string, a gene set name>: {
                "geneset_name": <string, a gene set name>,
                "geneset_description": <a string or None>,
                "genes": [
                    {
                        "gene_symbol": <string, a gene symbol or name>,
                        "gene_description": <a string or None>
                    },
                    ...
                ]
            },
            ...
        }
    """

    class myDialect(csv.excel):
        skipinitialspace = True

    def just(n, seq):
        it = iter(seq)
        for _ in range(n - 1):
            yield next(it, "")
        yield tuple(it)

    messagefn = context["messagefn"] if context else (lambda x: None)

    reader = csv.reader(f, dialect=myDialect())
    genesets = {}
    haveReadHeader = False
    lineno = 0
    for row in reader:
        lineno += 1
        # ignore empty rows
        if len(row) == 0:
            continue
        # if row starts with '#' it is a comment
        if row[0].startswith("#"):
            continue
        # if this is the first non-comment row, assume it is a header
        if not haveReadHeader:
            if row != Annotations.Genesets_Header:
                raise AnnotationsError("Geneset CSV file missing the required column header.")
            haveReadHeader = True
            continue

        geneset_name, geneset_description, gene_symbol, gene_description, _ = just(5, row)
        if not geneset_name:
            raise AnnotationsError(f"Geneset CSV missing required geneset or gene name on line {lineno}")
        if (not gene_symbol) and gene_description:
            messagefn(f"Warning: Missing gene name in geneset name {geneset_name} on line {lineno}.")

        if geneset_name in genesets:
            gs = genesets[geneset_name]
        else:
            gs = genesets[geneset_name] = {
                "geneset_name": geneset_name,
                "geneset_description": geneset_description,
                "genes": [],
            }
        # Use first geneset_description with a value
        if not gs["geneset_description"] and geneset_description:
            gs["geneset_description"] = geneset_description
        # add the gene if the gene_symbol is defined
        if gene_symbol:
            gs["genes"].append(
                {
                    "gene_symbol": gene_symbol,
                    "gene_description": gene_description,
                }
            )

    return genesets
