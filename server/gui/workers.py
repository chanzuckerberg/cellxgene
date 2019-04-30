import traceback

from PyQt5.QtCore import QRunnable, pyqtSlot

from server.gui.utils import WorkerSignals


class DataLoadWorker(QRunnable):
    def __init__(self, data_file, layout="umap", *args, **kwargs):
        super(DataLoadWorker, self).__init__()
        self.data_file = data_file
        self.layout = layout
        self.signals = WorkerSignals()

    @pyqtSlot()
    def run(self):
        if not self.data_file:
            self.signals.finished.emit()
            return

        # delayed import to speed load
        from server.app.scanpy_engine.scanpy_engine import ScanpyEngine
        args = {
            "layout": self.layout,
            "diffexp": "ttest",
            "max_category_items": 100,
            "diffexp_lfc_cutoff": 0.01,
            "obs_names": None,
            "var_names": None,
        }
        try:
            data_results = ScanpyEngine(self.data_file, args)
        except Exception as e:
            traceback.print_exc()
            self.signals.error.emit(str(e))
        else:
            self.signals.result.emit(data_results)
        finally:
            self.signals.finished.emit()


class ServerRunWorker(QRunnable):
    def __init__(self, app, host, port, *args, **kwargs):
        super(ServerRunWorker, self).__init__()
        self.app = app
        self.host = host
        self.port = port

    @pyqtSlot()
    def run(self):
        self.app.run(host=self.host, debug=False, port=self.port, threaded=True)
