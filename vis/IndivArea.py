from PySide2.QtCore import *
from PySide2.QtGui import *
from PySide2.QtWidgets import *

from Utils import Utils

class IndivArea(QWidget):
    def __init__(self, data_tools, initial_pop_index, *args, **kwargs):
        QWidget.__init__(self, *args, **kwargs)
        self.data_tools = data_tools
        self.pop_index = initial_pop_index

        export_button = QPushButton('Export')
        export_button.setSizePolicy(QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed))
        export_button.clicked.connect(self.export_gene_desc)

        data, max_cols = self.data_tools.get_gene_descs(self.pop_index)
        self.table = QTableWidget(len(data) - 1, max_cols + 1) #first row is headers, length is max_cols + 1 (extra 1 is for the ea_step)
        headers = data[0]
        #headers.extend([''] * (len(headers) - max_cols))
        self.table.setHorizontalHeaderLabels(headers)
        self.table.verticalHeader().setVisible(False) #don't show row numbers along the left side
        self.populate_table(data, max_cols)
        
        layout = QVBoxLayout()
        layout.addWidget(self.table)
        layout.addWidget(export_button)
        self.setLayout(layout)

    def populate_table(self, data, max_cols):
        widths = list(map(len, data[0])) #start with widths of headers (note: length of data[0] (headers) is max_cols + 1)
        for i in range(1, len(data)): #start at 1 to skip over header row
            for j in range(max_cols + 1):
                if j < len(data[i]):
                    item = QTableWidgetItem(data[i][j])
                    self.table.setItem(i - 1, j, item)
                else:
                    item = QTableWidgetItem('')
                    self.table.setItem(i - 1, j, item)

                widths[j] = max(widths[j], len(item.text()))

        for col in range(max_cols + 1):
            self.table.setColumnWidth(col, widths[col] * 7)
        
    @Slot()
    def update_pop_index(self, pop_index):
        self.pop_index = pop_index
        self.table.setRowCount(0) #this deletes all existing rows
        data, max_cols = self.data_tools.get_gene_descs(pop_index)
        self.populate_table(data, max_cols)

    @Slot()
    def export_gene_desc(self):
        filename = Utils.run_save_csv_dialog()
        if filename:
            self.data_tools.export_gene_descs(self.pop_index, filename)
