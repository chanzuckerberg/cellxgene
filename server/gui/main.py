# flake8: noqa F403, F405
import sys
from os.path import splitext, basename

from cefpython3 import cefpython as cef
from PyQt5.QtCore import *
from PyQt5.QtWidgets import *

from server.gui.browser import CefWidget, CefApplication
from server.gui.cxg_server import cellxgeneServer, DataLoadWorker, ServerRunWorker
from server.gui.utils import WINDOWS, LINUX, MAC

# Configuration
# TODO remember this or calculate it?
WIDTH = 1024
HEIGHT = 768

# noinspection PyUnresolvedReferences
class MainWindow(QMainWindow):
    def __init__(self):
        super(MainWindow, self).__init__(None)
        self.thread_pool = QThreadPool()
        self.cef_widget = None
        self.data_widget = None
        self.server = cellxgeneServer(self)
        self.server.setup_app()
        self.run_server()
        self.setWindowTitle("cellxgene")

        # Strong focus - accepts focus by tab & click
        self.setFocusPolicy(Qt.StrongFocus)
        self.setup_layout()

    def setup_layout(self):
        self.resize(WIDTH, HEIGHT)
        self.cef_widget = CefWidget(self)
        # TODO rename navigation bar
        self.data_widget = LoadWidget(self)
        layout = QGridLayout()
        layout.addWidget(self.cef_widget, 1, 0)
        layout.addWidget(self.data_widget, 0, 0)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        layout.setRowStretch(0, 0)
        layout.setRowStretch(1, 1)
        frame = QFrame()
        frame.setLayout(layout)
        self.setCentralWidget(frame)

        if WINDOWS:
            # On Windows with PyQt5 main window must be shown first
            # before CEF browser is embedded, otherwise window is
            # not resized and application hangs during resize.
            self.show()

        # Browser can be embedded only after layout was set up
        self.cef_widget.embedBrowser()

        if LINUX:
            # On Linux with PyQt5 the QX11EmbedContainer widget is
            # no more available. An equivalent in Qt5 is to create
            # a hidden window, embed CEF browser in it and then
            # create a container for that hidden window and replace
            # cef widget in the layout with the container.
            self.container = QWidget.createWindowContainer(
                self.cef_widget.hidden_window, parent=self)
            layout.addWidget(self.container, 1, 0)

    def closeEvent(self, event):
        # Close browser (force=True) and free CEF reference
        if self.cef_widget.browser:
            self.cef_widget.browser.CloseBrowser(True)
            self.clear_browser_references()

    def run_server(self):
        worker = ServerRunWorker(self.server.app, host=self.server.host, port=self.server.port)
        self.thread_pool.start(worker)

    def clear_browser_references(self):
        # Clear browser references that you keep anywhere in your
        # code. All references must be cleared for CEF to shutdown cleanly.
        self.cef_widget.browser = None


# TODO make central location for methods?
MODES = ["umap", "tsne", "draw_graph_fa", "draw_graph_fr", "diffmap", "phate"]
class LoadWidget(QFrame):
    def __init__(self, parent):
        super(LoadWidget, self).__init__(parent=parent)

        # Init layout
        load_ui_layout = QVBoxLayout()

        layout = QGridLayout()
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        self.title = ""
        self.label = QLabel("cellxgene")
        layout.addWidget(self.label, 0, 0, 1, 3)
        self.load = QPushButton("Open...")
        self.load.clicked.connect(self.on_load)
        layout.addWidget(self.load, 1, 2)
        self.embedding_label = QLabel("Embedding: ")
        layout.addWidget(self.embedding_label, 1, 0)
        self.embeddings = QComboBox(self)
        self.embeddings.currentIndexChanged.connect(self.update_embedding)
        self.embeddings.addItems(MODES)
        self.embedding_selection = MODES[0]
        layout.addWidget(self.embeddings, 1, 1)

        self.error_label = QLabel("")
        layout.addWidget(self.error_label, 2, 1)
        # Layout
        self.setLayout(layout)

    @pyqtSlot(int)
    def update_embedding(self, idx):
        self.embedding_selection = MODES[idx]

    @pyqtSlot()
    def on_load(self):
        options = QFileDialog.Options()
        # options |= QFileDialog.DontUseNativeDialog
        file_name, _ = QFileDialog.getOpenFileName(self,
                                                  "Open H5AD File", "", "H5AD Files (*.h5ad)", options=options)
        self.title = splitext(basename(file_name))[0]
        # self.cef_widget.parent.server.load_data(fileName)
        worker = DataLoadWorker(file_name, self.embedding_selection)
        worker.signals.result.connect(self.on_data_success)
        worker.signals.error.connect(self.on_data_error)
        self.window().thread_pool.start(worker)

    def on_data_success(self, data):
        self.window().server.attach_data(data, self.title)
        self.navigate_to_location()

    def on_data_error(self, err):
        self.error_label.setText(err)

    def navigate_to_location(self, location="http://localhost:8000/"):
        self.window().cef_widget.browser.Navigate(location)

# TODO make this cleaner
# Document more
# rename?
def main():
    sys.excepthook = cef.ExceptHook  # To shutdown all CEF processes on error
    settings = {}
    # Instead of timer loop
    if MAC:
        settings["external_message_pump"] = True

    cef.Initialize(settings)
    app = CefApplication(sys.argv)
    main_window = MainWindow()
    main_window.show()
    main_window.activateWindow()
    main_window.raise_()
    app.exec_()
    if not cef.GetAppSetting("external_message_pump"):
        app.stopTimer()
    del main_window  # Just to be safe, similarly to "del app"
    del app  # Must destroy app object before calling Shutdown
    cef.Shutdown()
    sys.exit(0)

if __name__ == '__main__':
    main()
