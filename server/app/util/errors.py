class FilterError(Exception):
    """
    Raised when filter is malformed
    """

    def __init__(self, message):
        self.message = message


class InteractiveError(Exception):
    """
    Raised when computation would exceed interactive time
    """

    def __init__(self, message):
        self.message = message


class MimeTypeError(Exception):
    """
    Raised when incompatible MIME type selected
    """

    def __init__(self, message):
        self.message = message


class PrepareError(Exception):
    """
    Raised when data is misprepared
    """

    def __init__(self, message):
        self.message = message


class ScanpyFileError(Exception):
    """
    Raised when file loaded into scanpy is misformatted
    """

    def __init__(self, message):
        self.message = message
