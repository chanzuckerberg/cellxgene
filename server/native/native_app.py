
import sys
import threading
from cefpython3 import cefpython as cef
import ctypes
import platform

from PyQt5.QtGui import *
from PyQt5.QtCore import *
from PyQt5.QtWidgets import *

WindowUtils = cef.WindowUtils()

# Platforms
WINDOWS = (platform.system() == "Windows")
LINUX = (platform.system() == "Linux")
MAC = (platform.system() == "Darwin")

# Configuration
WIDTH = 1024
HEIGHT = 768

# OS differences
CefWidgetParent = QWidget
if LINUX:
    CefWidgetParent = QX11EmbedContainer

def main():
    sys.excepthook = cef.ExceptHook  # To shutdown all CEF processes on error
    settings = {}
    if MAC:
        settings["external_message_pump"] = True

    datad = DataDaemon()
    datad.start()

    cef.Initialize(settings)
    app = CefApplication(sys.argv)
    main_window = MainWindow(datad)
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

class DataDaemon():
    def __init__(self):
        self.httpd = None
        self.app = None

    def data_server(self):
        host = "127.0.0.1"
        port = 8000
        debug = False
        # TODO check if thread should have access to self?
        self.app.run(host=host, debug=debug, port=port, threaded=True)

    def start(self):
        from server.app.app import Server
        server = Server()
        self.app = server.create_app()
        self.app.config.update(DATASET_TITLE="DEMO!")

        # start server
        self.httpd = threading.Thread(target=self.data_server, daemon=True)
        self.httpd.start()

    def load_data(self, data):
        print("DataDaemon::load ", data)
        from server.app.scanpy_engine.scanpy_engine import ScanpyEngine
        args = {
            "layout": "umap",
            "diffexp": "ttest",
            "max_category_items": 100,
            "diffexp_lfc_cutoff": 0.01,
            "obs_names": None,
            "var_names": None,
        }
        self.app.data = ScanpyEngine(data, args)


class MainWindow(QMainWindow):
    def __init__(self, datad):
        super(MainWindow, self).__init__(None)
        self.cef_widget = None
        self.navigation_bar = None
        # TODO what is datad?
        self.datad = datad
        self.setWindowTitle("example")
        # TODO check this
        self.setFocusPolicy(Qt.StrongFocus)
        self.setupLayout()

    def setupLayout(self):
        self.resize(WIDTH, HEIGHT)
        self.cef_widget = CefWidget(self)
        self.navigation_bar = LoadWidget(self.cef_widget)
        layout = QGridLayout()
        layout.addWidget(self.cef_widget, 1, 0)
        layout.addWidget(self.navigation_bar, 0, 0)
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
            # noinspection PyUnresolvedReferences, PyArgumentList
            self.container = QWidget.createWindowContainer(
                    self.cef_widget.hidden_window, parent=self)
            # noinspection PyArgumentList
            layout.addWidget(self.container, 1, 0)

    def closeEvent(self, event):
        # Close browser (force=True) and free CEF reference
        if self.cef_widget.browser:
            self.cef_widget.browser.CloseBrowser(True)
            self.clear_browser_references()

    def clear_browser_references(self):
        # Clear browser references that you keep anywhere in your
        # code. All references must be cleared for CEF to shutdown cleanly.
        self.cef_widget.browser = None


class CefWidget(CefWidgetParent):
    def __init__(self, parent=None):
        super(CefWidget, self).__init__(parent)
        self.parent = parent
        self.browser = None
        # TODO why??
        self.hidden_window = None  # Required for PyQt5 on Linux
        self.show()

    # TODO focus why?
    def focusInEvent(self, event):
        # This event seems to never get called on Linux, as CEF is
        # stealing all focus due to Issue #284.
        if self.browser:
            if WINDOWS:
                WindowUtils.OnSetFocus(self.getHandle(), 0, 0, 0)
            self.browser.SetFocus(True)

    def focusOutEvent(self, event):
        # This event seems to never get called on Linux, as CEF is
        # stealing all focus due to Issue #284.
        if self.browser:
            self.browser.SetFocus(False)

    # TODO when does this happen?
    def embedBrowser(self):
        if LINUX:
            # noinspection PyUnresolvedReferences
            self.hidden_window = QWindow()
        window_info = cef.WindowInfo()
        rect = [0, 0, self.width(), self.height()]
        window_info.SetAsChild(self.getHandle(), rect)
        # TODO better splash
        self.browser = cef.CreateBrowserSync(window_info,
                                             url="http://localhost:8000/splash")
        # TODO is this necessary?
        # self.browser.SetClientHandler(LoadHandler())
        self.browser.SetClientHandler(FocusHandler(self))

    # TODO is this needed?
    def getHandle(self):
        if self.hidden_window:
            # PyQt5 on Linux
            return int(self.hidden_window.winId())
        else:
            return int(self.winId())


    def moveEvent(self, _):
        self.x = 0
        self.y = 0
        if self.browser:
            if WINDOWS:
                WindowUtils.OnSize(self.getHandle(), 0, 0, 0)
            elif LINUX:
                self.browser.SetBounds(self.x, self.y,
                                       self.width(), self.height())
            self.browser.NotifyMoveOrResizeStarted()

    def resizeEvent(self, event):
        size = event.size()
        if self.browser:
            if WINDOWS:
                WindowUtils.OnSize(self.getHandle(), 0, 0, 0)
            elif LINUX:
                self.browser.SetBounds(self.x, self.y,
                                       size.width(), size.height())
            self.browser.NotifyMoveOrResizeStarted()


class CefApplication(QApplication):
    def __init__(self, args):
        super(CefApplication, self).__init__(args)
        # TODO is this the best way?
        # do we need to do this?
        if not cef.GetAppSetting("external_message_pump"):
            self.timer = self.createTimer()

    def createTimer(self):
        timer = QTimer()
        timer.timeout.connect(self.onTimer)
        timer.start(10)
        return timer

    def onTimer(self):
        cef.MessageLoopWork()

    def stopTimer(self):
        # Stop the timer after Qt's message loop has ended
        self.timer.stop()


class FocusHandler(object):
    # TODO why do we need this?
    def __init__(self, cef_widget):
        self.cef_widget = cef_widget

    def OnSetFocus(self, **_):
        pass

    def OnGotFocus(self, browser, **_):
        # Temporary fix no. 1 for focus issues on Linux (Issue #284)
        if LINUX:
            print("[qt.py] FocusHandler.OnGotFocus:"
                  " keyboard focus fix no. 1 (Issue #284)")
            browser.SetFocus(True)

class LoadWidget(QFrame):
    def __init__(self, cef_widget):
        super(LoadWidget, self).__init__()
        self.cef_widget = cef_widget

        # Init layout
        layout = QGridLayout()
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)

        self.load = QPushButton("load")
        self.load.clicked.connect(self.onLoad)
        layout.addWidget(self.load, 0, 0)

        # Layout
        self.setLayout(layout)

    def onLoad(self):
        print("load!")
        options = QFileDialog.Options()
        # options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getOpenFileName(self,
            "Open H5AD File", "","H5AD Files (*.h5ad)", options=options)
        if fileName:
            # TODO handle this better
            # TODO thread this
            self.cef_widget.parent.datad.load_data(fileName)
            self.load.setEnabled(False)
            # TODO instead create cef_widget
            self.cef_widget.browser.Navigate("http://localhost:8000/")

    def createButton(self, name):
        return QPushButton(name)

if __name__ == '__main__':
    main()
