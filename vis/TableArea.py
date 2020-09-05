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
        self.cur_index = initial_index
        data = self.data_tools.get_protein_info_for_indiv(initial_index, self.data_tools.tag_type_dict["IndivStateAfterBind"])
        
        vlayout = QVBoxLayout()

        headers = [
            'Type',
            'Tag',
            'Action',
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

        radio_layout = QVBoxLayout()
        self.after_bind_button = QRadioButton("After Bind", self)
        self.after_bind_button.setChecked(True)
        self.after_prod_button = QRadioButton("After Prod", self)
        radio_layout.addWidget(self.after_bind_button)
        radio_layout.addWidget(self.after_prod_button)
        self.after_bind_button.toggled.connect(lambda checked: self.refresh(self.cur_index) if checked else None)
        self.after_prod_button.toggled.connect(lambda checked: self.refresh(self.cur_index) if checked else None)
        
        search_hlayout = QHBoxLayout()
        search_label = QLabel('Search:')
        self.searchEntry = QLineEdit()
        self.searchEntry.textChanged.connect(self.proxyModel.setFilterWildcard)
        search_hlayout.addWidget(search_label)
        search_hlayout.addWidget(self.searchEntry)

        button_hlayout = QHBoxLayout()
        self.sel_all_button = QPushButton("Select All")
        self.desel_all_button = QPushButton("Deselect All")
        self.sel_all_button.clicked.connect(self.handle_sel_all)
        self.desel_all_button.clicked.connect(self.handle_desel_all)
        button_hlayout.addWidget(self.sel_all_button)
        button_hlayout.addWidget(self.desel_all_button)

        vlayout.addWidget(QLabel("Checkpoint:"))
        vlayout.addLayout(radio_layout)
        vlayout.addLayout(search_hlayout)
        vlayout.addLayout(button_hlayout)
        vlayout.addWidget(self.table)
        self.setLayout(vlayout)
        
    def getCheckedInfo(self):
        return self.model.getCheckedInfo()

    def getVisibleRows(self):
        visible_rows = []
        for row in range(self.model.rowCount()):
            index = self.model.index(row, 0)
            if self.proxyModel.filterAcceptsRow(row, index):
                visible_rows.append(row)

        return visible_rows

    @Slot()
    def reset_colour_picker(self, new_ea_step):
        self.model.reset_colour_picker()

    @Slot()
    def refresh(self, index):
        if self.after_bind_button.isChecked():
            tag_type = self.data_tools.tag_type_dict['IndivStateAfterBind']
        else:
            tag_type = self.data_tools.tag_type_dict['IndivStateAfterProd']
            
        data = self.data_tools.get_protein_info_for_indiv(index, tag_type)
        self.model.refresh(data)
        self.cur_index = index

    @Slot()
    def handle_sel_all(self):
        row_indices = self.getVisibleRows()
        self.model.setCheckedRows(row_indices)

    @Slot()
    def handle_desel_all(self):
        self.model.setCheckedRows([])

    @Slot()
    def handle_rowChanged(self, checkedInfo):
        self.checksChanged.emit(checkedInfo)
