import traceback

from PyQt5.QtCore import QRunnable, pyqtSlot

from server.gui.utils import WorkerSignals


class cellxgeneServer():
    def __init__(self, parent, host="127.0.0.1", port=8000):
        self.app = None
        self.host = host
        self.port = port

    def setup_app(self):
        from server.app.app import Server
        server = Server()
        self.app = server.create_app()

    def attach_data(self, data, title="Demo"):
        self.app.config.update(DATASET_TITLE=title)
        self.app.data = data


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
