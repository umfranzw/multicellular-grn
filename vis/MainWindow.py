from PySide2.QtCore import *
from PySide2.QtGui import *
from PySide2.QtWidgets import *

from DataTools import DataTools
from TreeTools import TreeTools
from CellArea import CellArea
from TableArea import TableArea
from ToolbarArea import ToolbarArea
from GraphicsArea import GraphicsArea
from FitnessArea import FitnessArea
from NeighbourArea import NeighbourArea
from IndivArea import IndivArea
from SettingsArea import SettingsArea, Settings

class MainWindow(QMainWindow):
    def __init__(self):
        QMainWindow.__init__(self)
        self.setWindowTitle("Vis")

        # Data
        self.data_tools = DataTools('data')
        self.tree_tools = TreeTools(self.data_tools)
        self.run = self.data_tools.get_run()
        
        # Tool Bar
        self.toolbar = ToolbarArea(self.run, self.data_tools)
        self.addToolBar(self.toolbar)

        # Table Area
        self.table_area = TableArea(self.data_tools, self.toolbar.getIndex())

        # cell info area
        self.cell_area = CellArea(self.data_tools)

        # graphics area
        self.graphics_area = GraphicsArea(self.data_tools, self.tree_tools, self.toolbar.getIndex())
        self.graphics_area.setSizePolicy(QSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding))

        #fitness area
        self.fitness_area = FitnessArea(self.data_tools)

        #neighbour interaction area
        self.neighbour_area = NeighbourArea(self.data_tools, self.toolbar.getIndex())

        self.indiv_area = IndivArea(self.data_tools, self.toolbar.getIndex()[1])
        
        #settings area
        self.settings_area = SettingsArea()
        
        self.toolbar.indexChanged.connect(self.table_area.refresh)
        self.toolbar.indexChanged.connect(self.refresh_graphics_area)
        self.toolbar.indexChanged.connect(lambda index: self.neighbour_area.update(index))
        self.toolbar.eaStepChanged.connect(self.indiv_area.update_ea_step)
        self.toolbar.showBestChanged.connect(self.show_best)
        self.toolbar.showBestChanged.connect(lambda checked: self.table_area.refresh(self.toolbar.getIndex()))
        self.toolbar.showBestChanged.connect(lambda checked: self.refresh_graphics_area(self.toolbar.getIndex()))
        self.toolbar.eaStepChanged.connect(self.table_area.reset_colour_picker)
        self.table_area.checksChanged.connect(self.refresh_graphics_area)
        self.graphics_area.selectionChanged.connect(self.relay_selection_changed)
        self.settings_area.settingChanged.connect(lambda name: self.refresh_graphics_area(self.toolbar.getIndex()))

        #combine everything into a series of tabs
        tabs = QTabWidget()
        tabs.addTab(self.table_area, "Protein Info")
        tabs.addTab(self.cell_area, "Cell Info")
        tabs.addTab(self.fitness_area, "Fitness Info")
        tabs.addTab(self.neighbour_area, "Neighbour Info")
        tabs.addTab(self.indiv_area, "Individual Info")
        tabs.addTab(self.settings_area, "Settings")
        
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

    #arg can be the current index (tuple) or checkedInfo (list)
    @Slot()
    def refresh_graphics_area(self, arg):
        if type(arg) == tuple:
            index = arg
            checkedInfo = self.table_area.getCheckedInfo()
            
        elif type(arg) == list:
            index = self.toolbar.getIndex()
            checkedInfo = arg

        self.graphics_area.refresh(index, checkedInfo)

    @Slot()
    def show_best(self, is_checked):
        self.blockSignals(True)
        
        if is_checked:
            self.prev_index = self.toolbar.getIndex()
            best_index = self.data_tools.get_run_best_index()
            self.toolbar.setIndex(best_index)
            self.toolbar.disable_spin_buttons()

        else:
            self.toolbar.enable_spin_buttons()
            self.toolbar.setIndex(self.prev_index)

        self.blockSignals(False)
