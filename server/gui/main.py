# flake8: noqa F403, F405
from os.path import splitext, basename
import sys
import threading

from cefpython3 import cefpython as cef
from PySide2.QtCore import *
from PySide2.QtWidgets import *

from server.app.app import Server
from server.gui.browser import CefWidget, CefApplication
from server.gui.workers import DataLoadWorker, ServerRunWorker
from server.gui.utils import WINDOWS, LINUX, MAC, FileLoadSignals
from server.utils.constants import MODES


# Configuration
# TODO remember this or calculate it?
WIDTH = 1024
HEIGHT = 768

# noinspection PyUnresolvedReferences
class MainWindow(QMainWindow):
    def __init__(self):
        super(MainWindow, self).__init__(None)
        self.cef_widget = None
        self.data_widget = None
        self.server = Server()
        self.server.create_app()
        self.runServer()
        self.setWindowTitle("cellxgene")

        # Strong focus - accepts focus by tab & click
        self.setFocusPolicy(Qt.StrongFocus)
        self.setupLayout()
        self.setupMenu()

    def setupLayout(self):
        self.resize(WIDTH, HEIGHT)
        self.cef_widget = CefWidget(self)
        self.data_widget = LoadWidget(self)
        self.stacked_layout = QStackedLayout()
        self.stacked_layout.addWidget(self.data_widget)
        self.stacked_layout.addWidget(self.cef_widget)
        main_layout = QVBoxLayout()
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.setSpacing(0)
        main_layout.addLayout(self.stacked_layout)
        frame = QFrame()
        frame.setLayout(main_layout)
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
            # no longer available. An equivalent in Qt5 is to create
            # a hidden window, embed CEF browser in it and then
            # create a container for that hidden window and replace
            # cef widget in the layout with the container.
            self.container = QWidget.createWindowContainer(
                self.cef_widget.hidden_window, parent=self)
            stacked_layout.addWidget(self.container, 1, 0)

    def setupMenu(self):
        main_menu = self.menuBar()
        file_menu = main_menu.addMenu('File')
        load_action = QAction("Load file...", self)
        load_action.setStatusTip("Load file")
        load_action.setShortcut("Ctrl+O")
        load_action.triggered.connect(self.showLoad)
        file_menu.addAction(load_action)

    def showLoad(self):
        self.stacked_layout.setCurrentIndex(0)

    def closeEvent(self, event):
        # Close browser (force=True) and free CEF reference
        if self.cef_widget.browser:
            self.cef_widget.browser.CloseBrowser(True)
            self.clearBrowserReferences()

    def runServer(self):
        worker = ServerRunWorker(self.server.app, host="127.0.0.1", port=8000)
        self.httpd = threading.Thread(target=worker.run, daemon=True)
        self.httpd.start()

    def clearBrowserReferences(self):
        # Clear browser references that you keep anywhere in your
        # code. All references must be cleared for CEF to shutdown cleanly.
        self.cef_widget.browser = None


class LoadWidget(QFrame):
    def __init__(self, parent):
        super(LoadWidget, self).__init__(parent=parent)
        # Init layout
        self.MAX_CONTENT_WIDTH = 500
        load_ui_layout = QVBoxLayout()
        h_margin = (WIDTH - self.MAX_CONTENT_WIDTH) // 2
        if h_margin < 10:
            h_margin = 10
        load_ui_layout.setContentsMargins(h_margin, 20, h_margin, 20)
        logo_layout = QHBoxLayout()
        logo_layout.setContentsMargins(0, 0, 0, 20)

        load_layout = QGridLayout()
        load_layout.setContentsMargins(0, 0, 0, 0)
        load_layout.setSpacing(0)
        message_layout = QHBoxLayout()
        message_layout.setContentsMargins(0, 0, 0, 0)

        self.title = ""
        self.label = QLabel("cellxgene")
        logo_layout.addWidget(self.label)

        # UI section
        # TODO add load spinner
        # TODO add cancel button to send back to browser (if available)
        self.embedding_label = QLabel("embedding: ")
        load_layout.addWidget(self.embedding_label, 0, 0)
        self.file_label = QLabel("file: ")
        load_layout.addWidget(self.file_label, 0, 1)
        self.embeddings = QComboBox(self)
        self.embeddings.currentIndexChanged.connect(self.updateEmbedding)
        self.embeddings.addItems(MODES)
        self.embedding_selection = MODES[0]
        load_layout.addWidget(self.embeddings, 1, 0)

        self.load = QPushButton("Open...")
        self.load.clicked.connect(self.onLoad)
        load_layout.addWidget(self.load, 1, 1)

        # Error section
        self.error_label = QLabel("")
        self.error_label.setWordWrap(True)
        self.error_label.setFixedWidth(self.MAX_CONTENT_WIDTH)
        message_layout.addWidget(self.error_label, alignment=Qt.AlignTop)

        # Layout
        for l in [logo_layout, load_layout, message_layout ]:
            load_ui_layout.addLayout(l)

        load_ui_layout.setStretch(2, 10)
        self.setLayout(load_ui_layout)

        self.signals = FileLoadSignals()
        self.signals.selectedFile.connect(self.createScanpyEngine)

    def updateEmbedding(self, idx):
        self.embedding_selection = MODES[idx]

    def createScanpyEngine(self, file_name):
        worker = DataLoadWorker(file_name, self.embedding_selection)
        worker.signals.result.connect(self.onDataSuccess)
        worker.signals.error.connect(self.onDataError)
        self.load_worker = threading.Thread(target=worker.run, daemon=True)
        self.load_worker.start()

    def onLoad(self):
        options = QFileDialog.Options()
        # options |= QFileDialog.DontUseNativeDialog
        file_name, _ = QFileDialog.getOpenFileName(self,
                                                  "Open H5AD File", "", "H5AD Files (*.h5ad)", options=options)
        self.title = splitext(basename(file_name))[0]
        if file_name:
            self.signals.selectedFile.emit(file_name)


    def onDataSuccess(self, data):
        self.window().server.attach_data(data, self.title)
        self.navigateToLocation()
        # Reveal browser
        self.window().stacked_layout.setCurrentIndex(1)

    def onDataError(self, err):
        self.error_label.setText(f"Error: {err}")
        self.error_label.resize(self.MAX_CONTENT_WIDTH, self.error_label.height())

    def navigateToLocation(self, location="http://localhost:8000/"):
        self.window().cef_widget.browser.Navigate(location)


def main():
    # This generates an error.log file on error
    sys.excepthook = cef.ExceptHook  # To shutdown all CEF processes on error
    settings = {}
    # Instead of timer loop
    if MAC:
        settings["external_message_pump"] = True

    # Create and launch cef browser and qt window
    cef.Initialize(settings)
    app = CefApplication(sys.argv)
    main_window = MainWindow()
    main_window.show()
    main_window.activateWindow()
    main_window.raise_()
    app.exec_()

    # Clean up on close
    if not cef.GetAppSetting("external_message_pump"):
        app.stopTimer()
    # TODO clean up threads when we switch threading model
    del main_window  # Just to be safe, similarly to "del app"
    del app  # Must destroy app object before calling Shutdown
    cef.Shutdown()
    sys.exit(0)


if __name__ == '__main__':
    main()
