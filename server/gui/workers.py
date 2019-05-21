from multiprocessing import Process
import traceback

s
class EmittingProcess(Process):
    def __init__(self, transport, *arg, **kwargs):
        super(EmittingProcess, self).__init__()
        self.transport = transport

    def emit(self, signal_name, *args):
        message = (signal_name, *args)
        self.transport.send(message)


class Worker(EmittingProcess):
    def __init__(self, transport, data_file, title, host, port, layout=["umap"], *args, **kwargs):
        super(Worker, self).__init__()
        self.data_file = data_file
        self.layout = layout
        self.title = title
        self.host = host
        self.port = port

    def run(self):
        from server.app.app import Server
        from server.app.scanpy_engine.scanpy_engine import ScanpyEngine
        server = Server()
        app = server.create_app()
        args = {
            "layout": self.layout,
            "diffexp": "ttest",
            "max_category_items": 100,
            "diffexp_lfc_cutoff": 0.01,
            "obs_names": None,
            "var_names": None,
        }
        data = ScanpyEngine(self.data_file, args)
        server.attach_data(data, self.title)
        self.emit("ready")
        # TODO listen for finished to kill?
        server.app.run(host=self.host, debug=True, port=self.port, threaded=True)


class ServerRunWorker():
    def __init__(self, app, host, port, *args, **kwargs):
        super(ServerRunWorker, self).__init__()
        self.app = app
        self.host = host
        self.port = port

    def run(self):
        self.app.run(host=self.host, debug=False, port=self.port, threaded=True)
