from PySide2.QtCore import Slot
from PySide2.QtGui import QKeySequence
from PySide2.QtWidgets import QMainWindow, QAction, QSpinBox, QToolBar, QLabel
from PySide2.QtGui import QPixmap, QImage

from DataTools import DataTools
from TreeTools import TreeTools

class MainWindow(QMainWindow):
    def __init__(self):
        QMainWindow.__init__(self)
        self.setWindowTitle("Vis")

        # Menu
        self.menu = self.menuBar()
        self.file_menu = self.menu.addMenu("File")

        # Exit QAction
        exit_action = QAction("Exit", self)
        exit_action.setShortcut(QKeySequence.Quit)
        exit_action.triggered.connect(self.close)

        self.file_menu.addAction(exit_action)

        # Data
        self.data_tools = DataTools('data')
        self.tree_tools = TreeTools(self.data_tools)
        self.run = self.data_tools.get_run()
        
        # Tool Bar
        self.toolbar = QToolBar()
        
        self.indivSpin = QSpinBox()
        self.indivSpin.setRange(1, self.run.pop_size)
        
        self.eaStepSpin = QSpinBox()
        self.eaStepSpin.setRange(0, self.run.ea_steps)
        #self.eaStepSpin.setSingleStep(self.run.step_interval.step)
        
        self.regStepSpin = QSpinBox()
        self.regStepSpin.setRange(1, self.run.reg_steps + 1)

        self.toolbar.addWidget(QLabel("Indiv:"))
        self.toolbar.addWidget(self.indivSpin)
        self.toolbar.addWidget(QLabel("EA Step:"))
        self.toolbar.addWidget(self.eaStepSpin)
        self.toolbar.addWidget(QLabel("Reg Step:"))
        self.toolbar.addWidget(self.regStepSpin)
        
        self.addToolBar(self.toolbar)

        # Image
        self.image_label = QLabel()
        
        self.setCentralWidget(self.image_label)

        # Window dimensions
        geometry = qApp.desktop().availableGeometry(self)
        self.setFixedSize(geometry.width() * 0.8, geometry.height() * 0.7)

        self.update_image()

    def update_image(self):
        index = (self.indivSpin.value(), self.eaStepSpin.value(), self.regStepSpin.value())
        raw = self.tree_tools.gen_image(self.data_tools, index)
        image = QImage.fromData(raw, format='png')
        pixmap = QPixmap.fromImage(image)
        self.image_label.setPixmap(pixmap)

