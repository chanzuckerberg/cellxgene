import os
import tempfile
import fsspec
from datetime import datetime


class DataLocator():
    """
    DataLocator is a simple wrapper around fsspec functionality, and provides a
    set of functions to encapsulate a data location (URI or path), interogate
    metadata about the object at that location (size, existance, etc) and
    access the underlying data.

    https://filesystem-spec.readthedocs.io/en/latest/index.html

    Example:
        dl = DataLocator("/tmp/foo.h5ad")
        if dl.exists():
            print(dl.size())
            with dl.open() as f:
                thecontents = f.read()

    DataLocator will accept a URI or native path.  Error handling is as defined
    in fsspec.

    """

    def __init__(self, uri_or_path):
        self.uri_or_path = uri_or_path
        self.protocol, self.path = DataLocator._get_protocol_and_path(uri_or_path)
        # work-around for LocalFileSystem not treating file: and None as the same scheme/protocol
        self.cname = self.path if self.protocol == 'file' else self.uri_or_path
        # will throw RuntimeError if the protocol is unsupported
        self.fs = fsspec.filesystem(self.protocol)

    @staticmethod
    def _get_protocol_and_path(uri_or_path):
        if "://" in uri_or_path:
            protocol, path = uri_or_path.split("://", 1)
            # windows!!!  Ignore single letter drive identifiers,
            # eg, G:\foo.txt
            if len(protocol) > 1:
                return protocol, path
        return None, uri_or_path

    def exists(self):
        return self.fs.exists(self.cname)

    def size(self):
        return self.fs.size(self.cname)

    def lastmodtime(self):
        """ return datetime object representing last modification time, or None if unavailable """
        info = self.fs.info(self.cname)
        if self.islocal() and info is not None:
            return datetime.fromtimestamp(info['mtime'])
        else:
            return getattr(info, 'LastModified', None)

    def abspath(self):
        """
        return the absolute path for the locator - only really does something
        for file: protocol, as all others are already absolute
        """
        if self.islocal():
            return os.path.abspath(self.path)
        else:
            return self.uri_or_path

    def isfile(self):
        return self.fs.isfile(self.cname)

    def open(self, *args):
        return self.fs.open(self.uri_or_path, *args)

    def islocal(self):
        return self.protocol is None or self.protocol == 'file'

    def local_handle(self):
        if self.islocal():
            return LocalFilePath(self.path)

        # if not local, create a tmp file system object to contain the data,
        # and clean it up when done.
        with self.open() as src, tempfile.NamedTemporaryFile(prefix="cellxgene_", delete=False) as tmp:
            tmp.write(src.read())
            tmp.close()
            src.close()
            tmp_path = tmp.name
            return LocalFilePath(tmp_path, delete=True)


class LocalFilePath():
    def __init__(self, tmp_path, delete=False):
        self.tmp_path = tmp_path
        self.delete = delete

    def __enter__(self):
        return self.tmp_path

    def __exit__(self, *args):
        if self.delete:
            os.unlink(self.tmp_path)
