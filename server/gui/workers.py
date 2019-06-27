from multiprocessing import Process
import time

import requests

from server.gui.utils import SiteReadySignals


class EmittingProcess(Process):
    def __init__(self, parent_conn, child_conn, *arg, **kwargs):
        super(EmittingProcess, self).__init__()
        self.parent_conn = parent_conn
        self.child_conn = child_conn

    def run(self):
        self.parent_conn.close()

    def emit(self, signal_name, *args):
        message = (signal_name, *args)
        self.child_conn.send(message)


class Worker(EmittingProcess):
    def __init__(self, parent_conn, child_conn, data_file, host, port, title, engine_options, *args, **kwargs):
        super(Worker, self).__init__(parent_conn, child_conn)
        self.data_file = data_file
        self.host = host
        self.port = port
        self.title = title
        self.engine_options = engine_options

    def run(self):
        super(Worker, self).run()
        if not self.data_file:
            self.emit("finished")
            return
        from server.app.app import Server
        from server.app.scanpy_engine.scanpy_engine import ScanpyEngine
        # create server
        try:
            server = Server()
            server.create_app()
        except Exception as e:
            self.emit("server_error", str(e))
            self.emit("finished")
            return
        # load data
        try:
            args = {
                "max_category_items": 100,
                "diffexp_lfc_cutoff": 0.01,
                "obs_names": None,
                "var_names": None,
            }
            args.update(self.engine_options)
            data = ScanpyEngine(self.data_file, args)
            server.attach_data(data, self.title)
            self.emit("ready")
        except Exception as e:
            self.emit("engine_error", str(e))
            self.emit("finished")
            return
        # launch server
        try:
            server.app.run(host=self.host, debug=False, port=self.port, threaded=True)
        except Exception as e:
            self.emit("server_error", str(e))
        finally:
            self.emit("finished")


class SiteReadyWorker:
    def __init__(self, location):
        super(SiteReadyWorker, self).__init__()
        self.signals = SiteReadySignals()
        self.location = location

    def run(self):
        session = requests.Session()
        for i in range(90):
            try:
                session.head(self.location)
                self.signals.ready.emit()
                break
            except requests.exceptions.ConnectionError:
                time.sleep(1)
            except Exception as e:
                self.signals.error.emit(str(e))
        self.signals.timeout.emit()
