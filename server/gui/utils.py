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
        ready
        error - `str` error message
        result - `object` data returned from processing, anything
    """
    finished = Signal()
    error = Signal(str)
    result = Signal(object)
    ready = Signal()


class SiteReadySignals(QObject):
    """
    Defines the signals available from a running worker thread.
    Supported signals are:
        timeout
        ready
        error - `str` error message
    """
    ready = Signal()
    timeout = Signal()
    error = Signal(str)


class FileLoadSignals(QObject):
    selectedFile = Signal(str)


class Emitter:
    def __init__(self, transport, signals):
        self.transport = transport
        self.signals = signals()

    def _emit(self, signature, args=None):
        if args:
            getattr(self.signals, signature).emit(args)
        else:
            getattr(self.signals, signature).emit()

    def run(self):
        while True:
            try:
                signature = self.transport.recv()
            except Exception:
                break
            else:
                self._emit(*signature)
