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
        self.graphics_area.selectionChanged.connect(self.relay_selection_changed)

        #combine everything into a series of tabs
        tabs = QTabWidget()
        tabs.addTab(self.table_area, "Protein Info")
        tabs.addTab(self.cell_area, "Cell Info")
        
        splitter = QSplitter()
        splitter.addWidget(self.graphics_area)
        splitter.addWidget(tabs)

        self.setCentralWidget(splitter)

    @Slot()
    def relay_selection_changed(self, cells):
        index = self.toolbar.getIndex()
        self.cell_area.refresh(cells, index)

    def sizeHint(self):
        return QSize(1200, 1000)

    def closeEvent(self, event):
        #note: need to disconnect this signal so it won't try to fire as the graphics area is being destroyed
        self.graphics_area.disconnect_signals()
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
