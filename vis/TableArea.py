from PySide2.QtCore import *
from PySide2.QtGui import *
from PySide2.QtWidgets import *
from TableModel import TableModel
from CustomSortFilterProxyModel import CustomSortFilterProxyModel

class TableArea(QWidget):
    checksChanged = Signal(list)
    
    def __init__(self, data_tools, initial_index, *args, **kwargs):
        QWidget.__init__(self, *args, **kwargs)
        self.data_tools = data_tools
        data = self.data_tools.get_protein_info_for_tree(initial_index)
        
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
        self.model.rowChanged.connect(self.handle_rowChanged)

        self.proxyModel = CustomSortFilterProxyModel(self)
        self.proxyModel.setSourceModel(self.model)
        self.table = QTableView()
        self.table.setModel(self.proxyModel)
        self.table.setSortingEnabled(True)

        hlayout = QHBoxLayout()
        search_label = QLabel('Search:')
        self.searchEntry = QLineEdit()
        self.searchEntry.textChanged.connect(self.proxyModel.setFilterWildcard)
        hlayout.addWidget(search_label)
        hlayout.addWidget(self.searchEntry)
        
        vlayout.addLayout(hlayout)
        vlayout.addWidget(self.table)
        self.setLayout(vlayout)
        
    def getCheckedInfo(self):
        return self.model.getCheckedInfo()

    @Slot()
    def refresh(self, index):
        data = self.data_tools.get_protein_info_for_tree(index)
        self.model.refresh(data)

    @Slot()
    def handle_rowChanged(self, checkedInfo):
        self.checksChanged.emit(checkedInfo)
