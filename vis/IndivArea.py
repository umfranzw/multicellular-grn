from PySide2.QtCore import *
from PySide2.QtGui import *
from PySide2.QtWidgets import *

from Utils import Utils

class IndivArea(QWidget):
    def __init__(self, data_tools, initial_ea_step, *args, **kwargs):
        QWidget.__init__(self, *args, **kwargs)
        self.data_tools = data_tools
        self.ea_step = initial_ea_step

        export_button = QPushButton('Export')
        export_button.setSizePolicy(QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed))
        export_button.clicked.connect(self.export_gene_desc)

        data, max_cols = self.data_tools.get_gene_descs(self.ea_step)
        self.table = QTableWidget(len(data) - 1, max_cols)
        headers = data[0]
        headers.extend([''] * (len(headers) - max_cols))
        self.table.setHorizontalHeaderLabels(headers)
        self.populate_table(data, max_cols)
        
        layout = QVBoxLayout()
        layout.addWidget(self.table)
        layout.addWidget(export_button)
        self.setLayout(layout)

    def populate_table(self, data, max_cols):
        widths = [0] * max_cols
        for i in range(1, len(data)): #start at 1 to skip over header row
            for j in range(max_cols):
                if j < len(data[i]):
                    item = QTableWidgetItem(data[i][j])
                    self.table.setItem(i, j, item)
                else:
                    item = QTableWidgetItem('')
                    self.table.setItem(i, j, item)

                widths[j] = max(widths[j], len(item.text()))

        for col in range(max_cols):
            self.table.setColumnWidth(col, widths[col] * 7)
        
    @Slot()
    def update_ea_step(self, ea_step):
        self.ea_step = ea_step
        self.table.setRowCount(0) #this deletes all existing rows
        data, max_cols = self.data_tools.get_gene_descs()
        self.populate_table(data, max_cols)

    @Slot()
    def export_gene_desc(self):
        filename = Utils.run_save_csv_dialog()
        if filename:
            self.data_tools.export_gene_descs(self.ea_step, filename)
