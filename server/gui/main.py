# flake8: noqa F403, F405
from functools import partialmethod
from multiprocessing import Pipe, Process, freeze_support
from os import environ
from os.path import splitext, basename, dirname, join
import sys
import threading

from cefpython3 import cefpython as cef
import PySide2
from PySide2.QtGui import *
from PySide2.QtCore import *
from PySide2.QtWidgets import *

import server.gui.cellxgene_rc
from server.gui.browser import CefWidget, CefApplication
from server.gui.workers import Worker, SiteReadyWorker
from server.gui.utils import WINDOWS, LINUX, MAC, FileLoadSignals, Emitter, WorkerSignals, FileChanged
from server.utils.utils import find_available_port

if WINDOWS or LINUX:
    dirname = dirname(PySide2.__file__)
    plugin_path = join(dirname, 'plugins', 'platforms')
    environ['QT_QPA_PLATFORM_PLUGIN_PATH'] = plugin_path

# Configuration
# TODO remember this or calculate it?
WIDTH = 1300
HEIGHT = 800
MAX_CONTENT_WIDTH = 700
GUI_PORT = find_available_port("localhost")
BROWSER_INDEX = 0
LOAD_INDEX = 1


class MainWindow(QMainWindow):
    def __init__(self):
        super(MainWindow, self).__init__(None)
        self.cef_widget = None
        self.data_widget = None
        self.stacked_layout = None
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
        self.clearMessages()
        self.stacked_layout.setCurrentIndex(LOAD_INDEX)

    def clearMessages(self):
        self.data_widget.reset()

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
        load_ui_layout = QVBoxLayout()
        h_margin = (WIDTH - MAX_CONTENT_WIDTH) // 2
        if h_margin < 10:
            h_margin = 10
        load_ui_layout.setContentsMargins(h_margin, 20, h_margin, 20)
        logo_layout = QHBoxLayout()
        logo_layout.setContentsMargins(0, 0, 0, 20)

        file_layout = QVBoxLayout()

        message_layout = QHBoxLayout()
        message_layout.setContentsMargins(0, 0, 0, 0)

        self.serverError = False
        self.file_name = FilePath()

        self.label = QLabel()
        logo = QPixmap(":/logo.png")
        self.label.setPixmap(logo)
        self.label.setContentsMargins(100, 0, 100, 0)
        logo_layout.addWidget(self.label)

        # UI section
        # TODO add cancel button to send back to browser (if available)

        self.file_area = FileArea(self)
        self.file_name.signals.changed.connect(self.updatePath)

        self.launch_widget = QLabel("Select a file to launch cellxgene")
        # self.launch_widget.setEnabled(False)
        # self.launch_widget.clicked.connect(self.onLoad)

        self.progress = QProgressBar()
        self.progress.setTextVisible(False)

        file_layout.addWidget(self.file_area)
        self.loading_layout = QStackedLayout()
        self.loading_layout.addWidget(self.launch_widget)
        self.loading_layout.addWidget(self.progress)
        file_layout.addLayout(self.loading_layout)
        file_layout.setStretch(0, 10)

        # Error section
        self.error_label = QLabel("")
        self.error_label.setWordWrap(True)
        self.error_label.setFixedWidth(MAX_CONTENT_WIDTH)
        message_layout.addWidget(self.error_label)

        # Options Form

        # Layout
        for l in [logo_layout, file_layout, message_layout]:
            load_ui_layout.addLayout(l)

        #TODO remove magic number
        load_ui_layout.setStretch(1, 10)
        self.setLayout(load_ui_layout)

        self.timer = QTimer()
        self.timer.setInterval(100)
        self.timer.timeout.connect(self.updateProgress)

        self.signals = FileLoadSignals()
        self.signals.selectedFile.connect(self.createScanpyEngine)
        self.signals.error.connect(self.onError)

    def updatePath(self):
        file_name = self.file_name.value
        if file_name:
            self.file_area.label.setText("File: " + file_name)
        else:
            self.file_area.label.setText("")
        # self.launch_widget.setEnabled(bool(file_name))

    def reset(self):
        self.loading_layout.setCurrentIndex(0)
        self.timer.stop()
        self.error_label.setText("")
        self.file_name.updateValue(None)

    def updateProgress(self):
        curr_val = self.progress.value()
        next_val = (curr_val + 1) % 100
        self.progress.setValue(next_val)

    def resetProgress(self):
        self.progress.setValue(0)
        self.loading_layout.setCurrentIndex(0)
        self.timer.stop()

    def createScanpyEngine(self, file_name):
        title = splitext(basename(file_name))[0]
        self.window().setupServer()
        worker = Worker(self.window().parent_conn, self.window().child_conn, file_name, host="127.0.0.1",
                        port=GUI_PORT, title=title, engine_options={})
        self.window().load_emitter.signals.ready.connect(self.onDataReady)
        self.window().load_emitter.signals.engine_error.connect(self.onServerError)
        self.window().load_emitter.signals.server_error.connect(self.onServerError)
        # Error is generic error from emitter
        self.window().load_emitter.signals.error.connect(self.onServerError)
        self.window().worker = Process(target=worker.run, daemon=True)
        self.window().worker.start()
        self.window().child_conn.close()

    def onLoad(self):
        if self.file_name.value:
            # Reset error on reload
            self.serverError = False
            self.loading_layout.setCurrentIndex(1)
            self.timer.start()
            self.signals.selectedFile.emit(self.file_name.value)
        else:
            self.signals.error.emit("Please select a file before launching.")

    def onDataReady(self):
        self.site_ready_worker = SiteReadyWorker(self.window().url)
        self.site_ready_worker.signals.ready.connect(self.onServerReady)
        self.site_ready_worker.signals.error.connect(self.onServerError)

        srw_thread = threading.Thread(target=self.site_ready_worker.run, daemon=True)
        srw_thread.start()

    def onServerReady(self):
        if not self.serverError:
            self.resetProgress()
            self.window().cef_widget.browser.Navigate(self.window().url)
            self.window().showBrowser()

    def onError(self, err, server_error=False):
        # Restart worker
        if server_error:
            self.serverError = True
            # Report error and switch to load screen
        self.window().shutdownServer()
        self.resetProgress()
        self.window().stacked_layout.setCurrentIndex(LOAD_INDEX)
        self.error_label.setText(f"Error: {err}")
        self.error_label.resize(MAX_CONTENT_WIDTH, self.error_label.height())
        self.window().repaint()

    onServerError = partialmethod(onError, server_error=True)

class FilePath(QObject):
    def __init__(self):
        super(FilePath, self).__init__()
        self.value = ""
        self.signals = FileChanged()

    def updateValue(self, path=None):
        self.value = path
        self.signals.changed.emit(self.value != path)


class FileArea(QFrame):
    def __init__(self, parent):
        super(FileArea, self).__init__()
        self.setFrameShape(QFrame.Box)
        self.setMinimumHeight(100)
        self.setFixedWidth(MAX_CONTENT_WIDTH)
        self.setAcceptDrops(True)
        self.instructions = QLabel(self)
        self.instructions.setText("Drag & Drop a h5ad file to load or open")
        self.instructions.setGeometry(10, 10, MAX_CONTENT_WIDTH, self.instructions.height())
        self.loadButton = QPushButton("Open...", parent=self)
        x_pos = (MAX_CONTENT_WIDTH - self.loadButton.width()) / 2
        self.loadButton.setGeometry(x_pos, 50, self.loadButton.width(), self.loadButton.height())
        self.loadButton.clicked.connect(self.fileBrowse)
        self.label = QLabel(self)
        self.label.setGeometry(10, 75, MAX_CONTENT_WIDTH, self.label.height())

    def fileBrowse(self):
        options = QFileDialog.Options()
        # options |= QFileDialog.DontUseNativeDialog
        file_name, _ = QFileDialog.getOpenFileName(self,
                                                   "Open H5AD File", "", "H5AD Files (*.h5ad)", options=options)
        if file_name:
            self.parent().file_name.updateValue(file_name)
            self.parent().onLoad()

    def dragEnterEvent(self, e):
        if e.mimeData().hasUrls:
            e.accept()
        else:
            e.ignore()

    def dragMoveEvent(self, e):
        if e.mimeData().hasUrls:
            e.accept()
        else:
            e.ignore()

    def dropEvent(self, e):
        """
        Drop files directly onto the widget
        File locations are stored in fname
        :param e:
        :return:
        """
        if e.mimeData().hasUrls:
            e.setDropAction(Qt.CopyAction)
            e.accept()
            for url in e.mimeData().urls():
                file_name = str(url.toLocalFile())
            self.parent().file_name.updateValue(file_name)
            self.parent().onLoad()
        else:
            e.ignore()


def main():
    freeze_support()
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
    main_window.setWindowTitle("cellxgene")
    main_window.setUnifiedTitleAndToolBarOnMac(True)
    main_window.setWindowIcon(QIcon(":icon.png"))
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
