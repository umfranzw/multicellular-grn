from PySide2.QtCore import *
from PySide2.QtGui import *
from PySide2.QtWidgets import *

from DataTools import DataTools
from TreeTools import TreeTools
from TableModel import TableModel
from CustomSortFilterProxyModel import CustomSortFilterProxyModel

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
        self.scene = QGraphicsScene()
        self.run = self.data_tools.get_run()
        
        # Tool Bar
        toolbar = self.build_toolbar()
        self.addToolBar(toolbar)

        # Table Area
        table_area = self.build_table()

        # Graphics View Area
        graphics_view = QGraphicsView(self.scene)
        #self.scene.selectionChanged.connect(lambda: self.update_cell_area(self.scene.selectedItems()))

        #combine the view and table areas horizontally
        centralWidget = QWidget(self)
        hlayout = QHBoxLayout()
        hlayout.addWidget(graphics_view)
        hlayout.addWidget(table_area)
        centralWidget.setLayout(hlayout)
        
        self.setCentralWidget(centralWidget)

        # Window dimensions
        geometry = qApp.desktop().availableGeometry(self)
        self.setFixedSize(geometry.width() * 0.8, geometry.height() * 0.7)

        self.update_view()

    def build_table(self):
        index = self.getIndex()
        data = self.data_tools.get_protein_info_for_tree(index)
        
        widget = QWidget(self)
        vlayout = QVBoxLayout()

        headers = [
            'Type',
            'Fcn',
            'Action',
            'Loc',
            'Arg',
            'isInitial',
            'Colour',
        ]
        self.model = TableModel(data, headers)
        self.model.rowChanged.connect(self.update_view)

        self.proxyModel = CustomSortFilterProxyModel(self)
        self.proxyModel.setSourceModel(self.model)
        self.table = QTableView()
        self.table.setModel(self.proxyModel)
        self.table.setSortingEnabled(True)

        self.searchEntry = QLineEdit()
        self.searchEntry.textChanged.connect(self.proxyModel.setFilterWildcard)
        
        vlayout.addWidget(self.searchEntry)
        vlayout.addWidget(self.table)
        widget.setLayout(vlayout)

        return widget

    def build_toolbar(self):
        toolbar = QToolBar()
        
        self.eaStepSpin = QSpinBox()
        self.eaStepSpin.setRange(0, self.run.ea_steps)
        self.eaStepSpin.valueChanged.connect(self.update_table)
        self.eaStepSpin.valueChanged.connect(self.update_view)
        #self.eaStepSpin.setSingleStep(self.run.step_interval.step)

        self.indivSpin = QSpinBox()
        self.indivSpin.setRange(1, self.run.pop_size)
        self.indivSpin.valueChanged.connect(self.update_table)
        self.indivSpin.valueChanged.connect(self.update_view)
        
        self.regStepSpin = QSpinBox()
        self.regStepSpin.setRange(1, self.run.reg_steps + 1)
        self.regStepSpin.valueChanged.connect(self.update_table)
        self.regStepSpin.valueChanged.connect(self.update_view)

        toolbar.addWidget(QLabel("EA Step:"))
        toolbar.addWidget(self.eaStepSpin)
        toolbar.addWidget(QLabel("Indiv:"))
        toolbar.addWidget(self.indivSpin)
        toolbar.addWidget(QLabel("Reg Step:"))
        toolbar.addWidget(self.regStepSpin)
        
        return toolbar

    def closeEvent(self, event):
        self.data_tools.close()

    def getIndex(self):
        return (self.eaStepSpin.value(), self.indivSpin.value(), self.regStepSpin.value())
    
    @Slot()
    def update_view(self):
        self.tree_tools.clear_scene(self.scene) #clear any previous tree
        index = self.getIndex()
        checked_info = self.model.getCheckedInfo() #[(julia_props_obj, QColor), ...]
        self.tree_tools.draw_scene(self.scene, index, checked_info)

    #note: there should only be at most one element in graph_items, since only one cell can be selected at once
    @Slot()
    def update_table(self):
        index = self.getIndex()
        data = self.data_tools.get_protein_info_for_tree(index)
        self.model.refresh(data)
