class FilterError(Exception):
    """
    Raised when filter is malformed
    """

    pass


class JSONEncodingValueError(Exception):
    """
    Raised when data cannot be encoded into json
    """

    pass


class MimeTypeError(Exception):
    """
    Raised when incompatible MIME type selected
    """

    pass


class PrepareError(Exception):
    """
    Raised when data is misprepared
    """

    pass


class DatasetAccessError(Exception):
    """
    Raised when file loaded into a DataAdaptor is misformatted
    """

    pass


class DisabledFeatureError(Exception):
    """
    Raised when an attempt to use a disabled feature occurs
    """

    pass


class AnnotationsError(Exception):
    """
    Raised when an attempt to use the annotations feature fails
    """

    pass


class OntologyLoadFailure(Exception):
    """
    Raised when reading the ontology file fails
    """

    pass


class ConfigurationError(Exception):
    """
    Raised when checking configuration errors
    """

    pass


class ExceedsLimitError(Exception):
    """
    Raised when an HTTP request exceeds a limit/quota
    """

    pass


class ComputeError(Exception):
    """
    Raised when an error occurs during a compute algorithm (such as diffexp)
    """

    pass
