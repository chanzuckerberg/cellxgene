import platform

from PyQt5.QtCore import QObject, pyqtSignal

WINDOWS = (platform.system() == "Windows")
LINUX = (platform.system() == "Linux")
MAC = (platform.system() == "Darwin")


class WorkerSignals(QObject):
    '''
    Defines the signals available from a running worker thread.
    Supported signals are:
    finished
    No data
    error
    `tuple` (exctype, value, traceback.format_exc() )
    result
    `object` data returned from processing, anything
    '''
    finished = pyqtSignal()
    error = pyqtSignal(tuple)
    result = pyqtSignal(object)
