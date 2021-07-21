from http import HTTPStatus


class CellxgeneException(Exception):
    """Base class for cellxgene exceptions"""

    def __init__(self, message):
        self.message = message
        super().__init__(message)


class RequestException(CellxgeneException):
    """Baseclass for exceptions that can be raised from a request."""

    # The default status code is 400 (Bad Request)
    default_status_code = HTTPStatus.BAD_REQUEST

    def __init__(self, message, status_code=None):
        super().__init__(message)
        self.status_code = status_code or self.default_status_code


def define_exception(name, doc):
    globals()[name] = type(name, (CellxgeneException,), dict(__doc__=doc))


def define_request_exception(name, doc, default_status_code=HTTPStatus.BAD_REQUEST):
    globals()[name] = type(name, (RequestException,), dict(__doc__=doc, default_status_code=default_status_code))


define_request_exception("FilterError", "Raised when filter is malformed")
define_request_exception("JSONEncodingValueError", "Raised when data cannot be encoded into json")
define_request_exception("MimeTypeError", "Raised when incompatible MIME type selected")
define_request_exception("DatasetAccessError", "Raised when file loaded into a DataAdaptor is misformatted")
define_request_exception("DisabledFeatureError", "Raised when an attempt to use a disabled feature occurs")
define_request_exception("AnnotationsError", "Raised when an attempt to use the annotations feature fails")
define_request_exception(
    "ComputeError",
    "Raised when an error occurs during a compute algorithm (such as diffexp)",
    HTTPStatus.INTERNAL_SERVER_ERROR,
)
define_request_exception("ExceedsLimitError", "Raised when an HTTP request exceeds a limit/quota")
define_request_exception("ColorFormatException", "Raised when color helper functions encounter an unknown color format")
define_request_exception(
    "AuthenticationError", "Raised when there is an authentication error", default_status_code=HTTPStatus.UNAUTHORIZED
)

define_request_exception(
    "AnnotationCategoryNameError",
    "Raised when an annotation category name cant be saved",
    default_status_code=HTTPStatus.UNPROCESSABLE_ENTITY,
)

define_exception("ConfigurationError", "Raised when checking configuration errors")
define_exception("PrepareError", "Raised when data is misprepared")
define_exception("SecretKeyRetrievalError", "Raised when get_secret_key from AWS fails")
define_exception("ObsoleteRequest", "Raised when the request is no longer valid.")
define_exception("UnsupportedSummaryMethod", "Raised when a gene set summary method is unknown or unsupported.")
define_request_exception("DatasetNotFoundException", "Raised when the dataset was not found in the data portal")
