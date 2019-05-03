import platform

from PySide2.QtCore import QObject, Signal

# Detect OS
WINDOWS = (platform.system() == "Windows")
LINUX = (platform.system() == "Linux")
MAC = (platform.system() == "Darwin")


class WorkerSignals(QObject):
    """
    Defines the signals available from a running worker thread.
    Supported signals are:
        finished
        error - `str` error message
        result - `object` data returned from processing, anything
    """
    finished = Signal()
    error = Signal(str)
    result = Signal(object)


class FileLoadSignals(QObject):
    selectedFile = Signal(str)
