from http import HTTPStatus


class RequestException(Exception):
    """Baseclass for exceptions that can be raised from a request."""

    # The default status code is 400 (Bad Request)
    default_status_code = HTTPStatus.BAD_REQUEST

    def __init__(self, message, status_code=None):
        Exception.__init__(self)
        self.message = message
        self.status_code = status_code or self.default_status_code


class FilterError(RequestException):
    """Raised when filter is malformed"""

    pass


class JSONEncodingValueError(RequestException):
    """Raised when data cannot be encoded into json"""

    pass


class MimeTypeError(RequestException):
    """Raised when incompatible MIME type selected"""

    pass


class DatasetAccessError(RequestException):
    """Raised when file loaded into a DataAdaptor is misformatted"""

    pass


class DisabledFeatureError(RequestException):
    """Raised when an attempt to use a disabled feature occurs"""

    pass


class AnnotationsError(RequestException):
    """Raised when an attempt to use the annotations feature fails"""

    pass


class ComputeError(RequestException):
    """Raised when an error occurs during a compute algorithm (such as diffexp)"""

    default_status_code = HTTPStatus.INTERNAL_SERVER_ERROR


class ExceedsLimitError(RequestException):
    """Raised when an HTTP request exceeds a limit/quota"""

    pass


class ColorFormatException(RequestException):
    """Raised when color helper functions encounter an unknown color format"""

    pass


class OntologyLoadFailure(Exception):
    """Raised when reading the ontology file fails"""

    pass


class ConfigurationError(Exception):
    """Raised when checking configuration errors"""

    pass


class PrepareError(Exception):
    """Raised when data is misprepared"""

    pass
