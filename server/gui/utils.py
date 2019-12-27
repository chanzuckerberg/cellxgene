import errno
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
    engine_error = Signal(str)
    server_error = Signal(str)
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
    error = Signal(str)


class FileChanged(QObject):
    changed = Signal(bool)


class Emitter:
    def __init__(self, transport, signals):
        self.transport = transport
        self.signals = signals()

    def _emit(self, signature, args=None):
        if args is None:
            getattr(self.signals, signature).emit()
        else:
            getattr(self.signals, signature).emit(args)

    def run(self):
        while True:
            try:
                signature = self.transport.recv()
            except EOFError:
                # Server done
                break
            except OSError as e:
                if e.errno == errno.EBADF:
                    break
                else:
                    self.signals.error.emit(str(e))
                    break
            except Exception as e:
                self.signals.error.emit(str(e))
                break
            else:
                self._emit(*signature)
