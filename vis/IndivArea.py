from PySide2.QtCore import *
from PySide2.QtGui import *
from PySide2.QtWidgets import *

from Utils import Utils

class IndivArea(QWidget):
    def __init__(self, data_tools, toolbar, *args, **kwargs):
        QWidget.__init__(self, *args, **kwargs)
        self.data_tools = data_tools
        self.toolbar = toolbar

        export_button = QPushButton('Export Gene Descriptions')
        export_button.setSizePolicy(QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed))
        export_button.clicked.connect(self.export_gene_desc)

        layout = QVBoxLayout()
        layout.addWidget(export_button)
        self.setLayout(layout)

    @Slot()
    def export_gene_desc(self):
        filename = Utils.run_save_csv_dialog()
        if filename:
            self.data_tools.export_gene_desc(self.toolbar.getIndex(), filename)
