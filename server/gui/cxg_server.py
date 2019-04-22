import traceback
import sys

from PyQt5.QtCore import QObject, QRunnable, pyqtSlot

from server.gui.utils import WorkerSignals


class cellxgeneServer(QObject):
    def __init__(self, parent, host="127.0.0.1", port=8000):
        super(cellxgeneServer, self).__init__(parent=parent)
        self.app = None
        self.host = host
        self.port = port

    def setup_app(self):
        from server.app.app import Server
        server = Server()
        self.app = server.create_app()
        self.app.config.update(DATASET_TITLE="DEMO!")
        # Technically a race condition
        self.run_server()

    def run_server(self):
        worker = ServerRunWorker(self.app, host=self.host, port=self.port)
        self.parent().thread_pool.start(worker)

    def attach_data(self, data):
        self.app.data = data


class DataLoadWorker(QRunnable):
    def __init__(self, data, *args, **kwargs):
        super(DataLoadWorker, self).__init__()
        self.data = data
        self.signals = WorkerSignals()

    @pyqtSlot()
    def run(self):
        from server.app.scanpy_engine.scanpy_engine import ScanpyEngine
        args = {
            "layout": "umap",
            "diffexp": "ttest",
            "max_category_items": 100,
            "diffexp_lfc_cutoff": 0.01,
            "obs_names": None,
            "var_names": None,
        }
        try:
            data_results = ScanpyEngine(self.data, args)
        except Exception as e:
            traceback.print_exc()
            exctype, value = sys.exc_info()[:2]
            self.signals.error.emit((exctype, value, traceback.format_exc()))
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
