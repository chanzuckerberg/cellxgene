class OptionsError(Exception):
    """
    Raised when user specified cli options specified fail parsing
    """

    def __init__(self, message):
        self.message = message
