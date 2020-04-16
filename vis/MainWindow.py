from PySide2.QtCore import *
from PySide2.QtGui import *
from PySide2.QtWidgets import *

from DataTools import DataTools
from TreeTools import TreeTools
from CellArea import CellArea
from TableArea import TableArea
from ToolbarArea import ToolbarArea
from GraphicsArea import GraphicsArea

class MainWindow(QMainWindow):
    def __init__(self):
        QMainWindow.__init__(self)
        self.setWindowTitle("Vis")

        # Data
        self.data_tools = DataTools('data')
        self.tree_tools = TreeTools(self.data_tools)
        self.run = self.data_tools.get_run()
        
        # Tool Bar
        self.toolbar = ToolbarArea(self.run)
        self.addToolBar(self.toolbar)

        # Table Area
        self.table_area = TableArea(self.data_tools, self.toolbar.getIndex())

        # cell info area
        self.cell_area = CellArea(self.data_tools)

        # graphics area
        self.graphics_area = GraphicsArea(self.tree_tools, self.toolbar.getIndex())
        self.graphics_area.setSizePolicy(QSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding))
        
        self.toolbar.indexChanged.connect(self.table_area.refresh)
        self.toolbar.indexChanged.connect(self.refresh_graphics_area)
        self.table_area.checksChanged.connect(self.refresh_graphics_area)
        self.graphics_area.selectionChanged.connect(self.cell_area.refresh)

        #combine the various areas in a grid
        centralWidget = QWidget(self)
        grid_layout = QGridLayout()
        grid_layout.addWidget(self.graphics_area, 0, 0, 2, 1)
        grid_layout.addWidget(self.table_area, 0, 1)
        grid_layout.addWidget(self.cell_area, 1, 1)
        centralWidget.setLayout(grid_layout)
        
        self.setCentralWidget(centralWidget)

    def sizeHint(self):
        return QSize(1200, 1000)

    def closeEvent(self, event):
        self.data_tools.close()

    @Slot()
    def refresh_graphics_area(self, arg):
        if type(arg) == tuple:
            index = arg
            checkedInfo = self.table_area.getCheckedInfo()
            
        elif type(arg) == list:
            index = self.toolbar.getIndex()
            checkedInfo = arg

        self.graphics_area.refresh(index, checkedInfo)
