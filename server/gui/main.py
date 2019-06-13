# flake8: noqa F403, F405
from functools import partialmethod
from multiprocessing import Pipe, Process
from os import environ
from os.path import splitext, basename, dirname, join
import sys
import threading

from cefpython3 import cefpython as cef
import PySide2
from PySide2.QtGui import *
from PySide2.QtCore import *
from PySide2.QtWidgets import *

from server.gui.browser import CefWidget, CefApplication
from server.gui.workers import Worker, SiteReadyWorker
from server.gui.utils import WINDOWS, LINUX, MAC, FileLoadSignals, Emitter, WorkerSignals
from server.utils.constants import MODES
from server.utils.utils import find_available_port

if WINDOWS or LINUX:
    dirname = dirname(PySide2.__file__)
    plugin_path = join(dirname, 'plugins', 'platforms')
    environ['QT_QPA_PLATFORM_PLUGIN_PATH'] = plugin_path

# Configuration
# TODO remember this or calculate it?
WIDTH = 1024
HEIGHT = 768
GUI_PORT = find_available_port("localhost")
BROWSER_INDEX = 0
LOAD_INDEX = 1


class MainWindow(QMainWindow):
    def __init__(self):
        super(MainWindow, self).__init__(None)
        self.cef_widget = None
        self.data_widget = None
        self.parent_conn, self.child_conn = None, None
        self.load_emitter = None
        self.emitter_thread = None
        self.worker = None
        self.url = f"http://localhost:{GUI_PORT}/"
        self.setWindowTitle("cellxgene")

        # Strong focus - accepts focus by tab & click
        self.setFocusPolicy(Qt.StrongFocus)
        self.setupLayout()
        self.setupMenu()

    def showBrowser(self):
        self.stacked_layout.setCurrentIndex(BROWSER_INDEX)

    def restartOnError(self):
        self.window().shutdownServer()
        # close emitter on error/finished
        self.parent_conn, self.child_conn = Pipe()
        self.load_emitter = Emitter(self.parent_conn, WorkerSignals)
        self.emitter_thread = threading.Thread(target=self.load_emitter.run, daemon=True)
        self.emitter_thread.start()
        # send to load with error message?

    def setupLayout(self):
        self.resize(WIDTH, HEIGHT)
        self.cef_widget = CefWidget(self)
        self.cef_widget.setSizePolicy(QSizePolicy(QSizePolicy.MinimumExpanding, QSizePolicy.MinimumExpanding))
        self.data_widget = LoadWidget(self)
        self.stacked_layout = QStackedLayout()
        self.stacked_layout.addWidget(self.cef_widget)
        self.stacked_layout.addWidget(self.data_widget)
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
            self.stacked_layout.replaceWidget(self.cef_widget, self.container)
        self.stacked_layout.setCurrentIndex(LOAD_INDEX)

    def setupServer(self):
        self.shutdownServer()
        # close emitter on error/finished
        self.parent_conn, self.child_conn = Pipe()
        self.load_emitter = Emitter(self.parent_conn, WorkerSignals)
        self.emitter_thread = threading.Thread(target=self.load_emitter.run, daemon=True)
        self.emitter_thread.start()
        # send to load with error message?

    def shutdownServer(self):
        if self.worker:
            self.worker.terminate()
        if self.parent_conn:
            self.parent_conn.close()

    def setupMenu(self):
        # TODO add communication to subprocess on reload
        main_menu = self.menuBar()
        file_menu = main_menu.addMenu('File')
        load_action = QAction("Load file...", self)
        load_action.setStatusTip("Load file")
        load_action.setShortcut("Ctrl+O")
        load_action.triggered.connect(self.showLoad)
        file_menu.addAction(load_action)

    def showLoad(self):
        self.stacked_layout.setCurrentIndex(LOAD_INDEX)

    def closeEvent(self, event):
        # Close browser (force=True) and free CEF reference
        if self.cef_widget.browser:
            self.cef_widget.browser.CloseBrowser(True)
            self.clearBrowserReferences()

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

        load_layout = QFormLayout()
        # load_layout.setContentsMargins(0, 0, 0, 0)
        # load_layout.setSpacing(0)
        message_layout = QHBoxLayout()
        message_layout.setContentsMargins(0, 0, 0, 0)

        self.serverError = False
        self.engine_options = {}

        self.label = QLabel("cellxgene")
        logo_layout.addWidget(self.label)
        logo = QPixmap()
        logo.load(":/../../cellxgene_logo.png")
        load_layout.addWidget(QLabel(logo))
        # UI section
        # TODO add load spinner
        # TODO add cancel button to send back to browser (if available)

        # Options
        self.title_widget = QLineEdit()
        self.title_widget.setAccessibleName("title")
        self.embedding_widget = QLineEdit()
        self.embedding_widget.setAccessibleName("layout")
        self.obs_widget = QLineEdit()
        self.obs_widget.setAccessibleName("obs_names")
        self.var_widget = QLineEdit()
        self.var_widget.setAccessibleName("var_names")
        self.max_category_items_widget = QSpinBox()
        self.max_category_items_widget.setAccessibleName("max_category_items")
        self.max_category_items_widget.setRange(1, 100000)
        self.max_category_items_widget.setValue(100)
        self.diff_exp_lfc_cutoff_widget = QDoubleSpinBox()
        self.diff_exp_lfc_cutoff_widget.setAccessibleName("diffexp_lfc_cutoff")
        self.diff_exp_lfc_cutoff_widget.setRange(0, 1)
        self.diff_exp_lfc_cutoff_widget.setValue(0.1)
        self.diff_exp_lfc_cutoff_widget.setDecimals(3)
        self.diff_exp_lfc_cutoff_widget.setStepType(QAbstractSpinBox.AdaptiveDecimalStepType)
        self.load_widget = QPushButton("Open...")
        self.load_widget.clicked.connect(self.onLoad)

        self.title_widget.textChanged.connect(self.updateOpt)
        self.embedding_widget.textChanged.connect(self.updateOpt)
        self.obs_widget.textChanged.connect(self.updateOpt)
        self.var_widget.textChanged.connect(self.updateOpt)
        self.max_category_items_widget.valueChanged.connect(self.updateOpt)
        self.diff_exp_lfc_cutoff_widget.valueChanged.connect(self.updateOpt)

        title_label = QLabel("title:")
        title_label.setToolTip("Title to display (if omitted will use file name)")
        embedding_label = QLabel("embedding:")
        embedding_label.setToolTip("Layout name, eg, 'umap'")
        obs_label = QLabel("obs names:")
        obs_label.setToolTip("Name of annotation field to use for observations")
        var_label = QLabel("var names:")
        var_label.setToolTip("Name of annotation to use for variables")
        max_category_items_label = QLabel("max categories:")
        max_category_items_label.setToolTip("Limits the number of categorical annotation items displayed")
        diff_exp_lfc_label = QLabel("diffexp cutoff:")
        diff_exp_lfc_label.setToolTip(
            "Relative expression cutoff used when selecting top N differentially expressed genes")
        load_label = QLabel("load file:")
        load_label.setToolTip("File to load (h5ad format)")

        load_layout.addRow(title_label, self.title_widget)
        load_layout.addRow(embedding_label, self.embedding_widget)
        load_layout.addRow(obs_label, self.obs_widget)
        load_layout.addRow(var_label, self.var_widget)
        load_layout.addRow(max_category_items_label, self.max_category_items_widget)
        load_layout.addRow(diff_exp_lfc_label, self.diff_exp_lfc_cutoff_widget)
        load_layout.addRow(load_label, self.load_widget)

        # Error section
        self.error_label = QLabel("")
        self.error_label.setWordWrap(True)
        self.error_label.setFixedWidth(self.MAX_CONTENT_WIDTH)
        message_layout.addWidget(self.error_label, alignment=Qt.AlignTop)

        # Layout
        for l in [logo_layout, load_layout, message_layout]:
            load_ui_layout.addLayout(l)

        load_ui_layout.setStretch(2, 10)
        self.setLayout(load_ui_layout)

        self.signals = FileLoadSignals()
        self.signals.selectedFile.connect(self.createScanpyEngine)



    def updateOpt(self, value):
        # get values
        caller = self.sender().accessibleName()
        val = None
        # Get value depending on element type
        if caller == "title":
            val = self.title_widget.text()
        elif caller == "layout":
            val = self.embedding_widget.text()
        elif caller == "obs_names":
            val = self.obs_widget.text()
        elif caller == "var_names":
            val = self.var_widget.text()
        elif caller == "max_category_items":
            val = self.max_category_items_widget.value()
        elif caller == "diffexp_lfc_cutoff":
            val = self.diff_exp_lfc_cutoff_widget.value()

        # Update options or delete if missing
        if val:
            self.engine_options[caller] = value
        elif caller in self.engine_options:
            del self.engine_options[caller]

    def createScanpyEngine(self, file_name):
        self.window().setupServer()
        worker = Worker(self.window().parent_conn, self.window().child_conn, file_name, host="127.0.0.1",
                        port=GUI_PORT, engine_options=self.engine_options)
        self.window().load_emitter.signals.ready.connect(self.onDataReady)
        self.window().load_emitter.signals.engine_error.connect(self.onServerError)
        self.window().load_emitter.signals.server_error.connect(self.onServerError)
        # Error is generic error from emitter
        self.window().load_emitter.signals.error.connect(self.onServerError)
        self.window().worker = Process(target=worker.run, daemon=True)
        self.window().worker.start()
        self.window().child_conn.close()

    def onLoad(self):
        options = QFileDialog.Options()
        # options |= QFileDialog.DontUseNativeDialog
        file_name, _ = QFileDialog.getOpenFileName(self,
                                                   "Open H5AD File", "", "H5AD Files (*.h5ad)", options=options)
        if "title" not in self.engine_options or not self.engine_options["title"]:
            self.engine_options["title"] = splitext(basename(file_name))[0]
        if file_name:
            self.signals.selectedFile.emit(file_name)
            # Reset error on reload
            self.serverError = False

    def onDataReady(self):
        self.site_ready_worker = SiteReadyWorker(self.window().url)
        self.site_ready_worker.signals.ready.connect(self.onServerReady)
        self.site_ready_worker.signals.error.connect(self.onServerError)

        srw_thread = threading.Thread(target=self.site_ready_worker.run, daemon=True)
        srw_thread.start()

    def onServerReady(self):
        if not self.serverError:
            self.window().cef_widget.browser.Navigate(self.window().url)
            self.window().showBrowser()

    def onError(self, err, server_error=False):
        # Restart worker
        if server_error:
            self.serverError = True
            # Report error and switch to load screen
            self.window().shutdownServer()
        self.window().stacked_layout.setCurrentIndex(LOAD_INDEX)
        self.error_label.setText(f"Error: {err}")
        self.error_label.resize(self.MAX_CONTENT_WIDTH, self.error_label.height())

    onServerError = partialmethod(onError, server_error=True)


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
    try:
        app.exec_()
    except Exception as e:
        raise
    finally:
        # Clean up on close
        if not cef.GetAppSetting("external_message_pump"):
            app.stopTimer()

        main_window.shutdownServer()
        del main_window  # Just to be safe, similarly to "del app"
        del app  # Must destroy app object before calling Shutdown
        cef.Shutdown()
        sys.exit(0)


if __name__ == '__main__':
    main()
