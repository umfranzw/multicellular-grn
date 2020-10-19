from PySide2.QtCore import *
from PySide2.QtGui import *
from PySide2.QtWidgets import *

from Utils import Utils

class IndivArea(QWidget):
    def __init__(self, data_tools, initial_pop_index, *args, **kwargs):
        QWidget.__init__(self, *args, **kwargs)
        self.data_tools = data_tools
        self.pop_index = initial_pop_index

        layout = QVBoxLayout()
        tabs = QTabWidget(self)
        genome_tab = self.build_genome_tab()
        reg_sim_tab = self.build_reg_sim_tab()
        fitness_info_tab = self.build_fitness_info_tab()

        tabs.addTab(genome_tab, "Genome")
        tabs.addTab(reg_sim_tab, "Reg Sim")
        tabs.addTab(fitness_info_tab, "Fitness Info")

        layout.addWidget(tabs)
        self.setLayout(layout)

    def build_fitness_info_tab(self):
        self.fitness_info_table, tab = self.build_table_tab(
            lambda: self.data_tools.get_fitness_info(self.pop_index),
            self.export_fitness_info
        )

        return tab
        
    def build_reg_sim_tab(self):
        self.reg_sim_table, tab = self.build_table_tab(
            lambda: self.data_tools.get_reg_sim_info(self.pop_index),
            self.export_reg_sim_info
        )

        return tab

    def build_genome_tab(self):
        self.genome_table, tab = self.build_table_tab(
            lambda: self.data_tools.get_gene_descs(self.pop_index),
            self.export_gene_desc
        )

        return tab
    
    def build_table_tab(self, data_fcn, export_fcn):
        tab = QWidget()
        layout = QVBoxLayout()
        
        export_button = QPushButton('Export')
        export_button.setSizePolicy(QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed))
        export_button.clicked.connect(export_fcn)

        data, max_cols = data_fcn()
        headers = data[0]
        table = QTableWidget(len(data) - 1, len(headers))
        table.setHorizontalHeaderLabels(headers)
        table.verticalHeader().setVisible(False) #don't show row numbers along the left side
        self.populate_table(data, table)
        
        layout.addWidget(table)
        layout.addWidget(export_button)
        tab.setLayout(layout)

        return table, tab

    def populate_table(self, data, table):
        widths = list(map(len, data[0])) #start with widths of headers
        num_cols = len(data[0]) #data[0] is headers
        
        for i in range(1, len(data)): #start at 1 to skip over header row
            for j in range(num_cols):
                if j < len(data[i]): #indiv genomes may grow over time...
                    item = QTableWidgetItem(data[i][j])
                    table.setItem(i - 1, j, item)
                else:
                    item = QTableWidgetItem('')
                    table.setItem(i - 1, j, item)

                widths[j] = max(widths[j], len(item.text()))

        for col in range(num_cols):
            table.setColumnWidth(col, widths[col] * 7)

    @Slot()
    def update_pop_index(self, pop_index):
        self.pop_index = pop_index
        self.update_table(
            self.genome_table,
            lambda: self.data_tools.get_gene_descs(pop_index)
        )

        self.update_table(
            self.reg_sim_table,
            lambda: self.data_tools.get_reg_sim_info(pop_index)
        )

        self.update_table(
            self.fitness_info_table,
            lambda: self.data_tools.get_fitness_info(pop_index)
        )

    def update_table(self, table, data_fcn):
        data, max_cols = data_fcn()
        table.clearContents() #clear old data
        #reset dimensions
        table.setRowCount(len(data) - 1)
        num_cols = len(data[0]) #length of header row
        table.setColumnCount(num_cols)
        self.populate_table(data, table)

    @Slot()
    def export_gene_desc(self):
        self.export_table(lambda: self.data_tools.export_gene_descs(self.pop_index, filename))

    @Slot()
    def export_reg_sim_info(self):
        self.export_table(lambda: self.data_tools.export_reg_sim_info(self.pop_index, filename))

    @Slot()
    def export_fitness_info(self):
        self.export_table(lambda: self.data_tools.export_fitness_info(self.pop_index, filename))

    def export_table(self, export_fcn):
        filename = Utils.run_save_csv_dialog()
        if filename:
            export_fcn()
