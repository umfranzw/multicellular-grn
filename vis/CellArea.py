from PySide2.QtCore import *
from PySide2.QtGui import *
from PySide2.QtWidgets import *
from PySide2.QtCharts import QtCharts

from Utils import Utils

class CellArea(QWidget):
    def __init__(self, data_tools, *args, **kwargs):
        QWidget.__init__(self, *args, **kwargs)
        self.data_tools = data_tools
        self.interaction_cell = None
        self.index = None
        
        layout = QVBoxLayout()
        tabs = QTabWidget(self)
        probs_tab = self.build_probs_tab()
        sensors_tab = self.build_sensors_tab()
        gene_states_tab = self.build_gene_states_tab()
        interaction_tab = self.build_interaction_tab()

        tabs.addTab(probs_tab, "Sym Probs")
        tabs.addTab(sensors_tab, "Sensors")
        tabs.addTab(gene_states_tab, "Gene States")
        tabs.addTab(interaction_tab, "Interaction")

        layout.addWidget(tabs)

        self.setLayout(layout)

    def build_gene_states_tab(self):
        tab = QWidget()
        layout = QVBoxLayout()
        
        #build table
        num_bind_sites = self.data_tools.get_run().bind_sites_per_gene
        headers = [
            'Gene Index',
            'Bind Logic',
        ]
        for i in range(num_bind_sites):
            headers.append('Bind Site {}'.format(i + 1))
            headers.append('Bound Protein {}'.format(i + 1))

        for i in range(num_bind_sites):
            headers.append('Prod Site {}'.format(i + 1))
            headers.append('Prod Rate {}'.format(i + 1))

        self.gs_table = QTableWidget(0, len(headers))
        self.gs_table.setHorizontalHeaderLabels(headers)

        #build gene scores chart
        self.gene_scores_chart = QtCharts.QChart()
        self.gene_scores_chart.setTitle('Gene Scores')
        self.gene_scores_view = QtCharts.QChartView(self.gene_scores_chart)
        self.gene_scores_view.setRenderHint(QPainter.RenderHint.Antialiasing)

        #save button
        save_button = QPushButton('Save')
        save_button.setSizePolicy(QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed))
        save_button.clicked.connect(lambda: Utils.save_chart_view(self.gene_scores_view))
        
        #put everything together
        layout.addWidget(self.gs_table)
        layout.addWidget(self.gene_scores_view)
        layout.addWidget(save_button)
        tab.setLayout(layout)
        
        return tab

    def build_interaction_tab(self):
        tab = QWidget()
        layout = QVBoxLayout()
        
        button = QPushButton('Generate')
        button.setSizePolicy(QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed))
        button.clicked.connect(self.update_interaction_img)
        layout.addWidget(button)
        
        scroll_area = QScrollArea()
        self.interaction_img = QLabel()
        scroll_area.setWidget(self.interaction_img)
        layout.addWidget(scroll_area)

        save_button = QPushButton('Save')
        save_button.setSizePolicy(QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed))
        save_button.clicked.connect(lambda: Utils.save_pixmap(self.interaction_img.pixmap()))
        layout.addWidget(save_button)

        tab.setLayout(layout)
        
        return tab

    def build_probs_tab(self):
        tab = QWidget()
        layout = QVBoxLayout()
        self.probs_chart = QtCharts.QChart()
        self.probs_chart.legend().setVisible(False)
        self.probs_chart.setTitle('Symbol Probs')

        chart_view = QtCharts.QChartView(self.probs_chart)
        chart_view.setRenderHint(QPainter.RenderHint.Antialiasing)
        self.probs_chart.resize(300, 300)

        save_button = QPushButton('Save')
        save_button.setSizePolicy(QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed))
        save_button.clicked.connect(lambda: Utils.save_chart_view(chart_view))
        
        layout.addWidget(chart_view)
        layout.addWidget(save_button)
        tab.setLayout(layout)

        scroll_area = QScrollArea()
        scroll_area.setWidget(tab)
        
        return scroll_area

    def clear_chart(self, chart):
        chart.removeAllSeries()
        for axis in chart.axes(Qt.Horizontal | Qt.Vertical):
            chart.removeAxis(axis)

    def build_sensors_tab(self):
        tab = QWidget()
        chart_widget = QWidget()
        vbox_layout = QVBoxLayout()
        grid_layout = QGridLayout()
        chart_size = (200, 200)

        self.top_chart = QtCharts.QChart()
        self.top_chart.legend().setVisible(False)
        self.top_view = QtCharts.QChartView(self.top_chart)
        self.top_view.setRenderHint(QPainter.RenderHint.Antialiasing)
        self.top_chart.resize(*chart_size)

        self.right_chart = QtCharts.QChart()
        self.right_chart.legend().setVisible(False)
        self.right_view = QtCharts.QChartView(self.right_chart)
        self.right_view.setRenderHint(QPainter.RenderHint.Antialiasing)
        self.right_chart.resize(*chart_size)

        self.bottom_chart = QtCharts.QChart()
        self.bottom_chart.legend().setVisible(False)
        self.bottom_view = QtCharts.QChartView(self.bottom_chart)
        self.bottom_view.setRenderHint(QPainter.RenderHint.Antialiasing)
        self.bottom_chart.resize(*chart_size)

        self.left_chart = QtCharts.QChart()
        self.left_chart.legend().setVisible(False)
        self.left_view = QtCharts.QChartView(self.left_chart)
        self.left_view.setRenderHint(QPainter.RenderHint.Antialiasing)
        self.left_chart.resize(*chart_size)

        grid_layout.addWidget(self.top_view, 0, 1)
        grid_layout.addWidget(self.right_view, 1, 2)
        grid_layout.addWidget(self.bottom_view, 2, 1)
        grid_layout.addWidget(self.left_view, 1, 0)

        save_button = QPushButton('Save')
        save_button.setSizePolicy(QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed))
        save_button.clicked.connect(lambda: self.save_sensor_charts(chart_size))
        
        chart_widget.setLayout(grid_layout)
        scroll_area = QScrollArea()
        scroll_area.setWidget(chart_widget)

        vbox_layout.addWidget(scroll_area)
        vbox_layout.addWidget(save_button)
        tab.setLayout(vbox_layout)

        return tab

    @Slot()
    def save_sensor_charts(self, chart_size):
        width, height = chart_size
        pixmap = QPixmap(3 * width, 3 * height)
        pixmap.fill(Qt.white)

        top_pixmap = QPixmap(width, height)
        self.top_view.render(top_pixmap)

        bottom_pixmap = QPixmap(width, height)
        self.bottom_view.render(bottom_pixmap)

        left_pixmap = QPixmap(width, height)
        self.left_view.render(left_pixmap)

        right_pixmap = QPixmap(width, height)
        self.right_view.render(right_pixmap)

        painter = QPainter(pixmap)
        painter.drawPixmap(QRect(width, 0, width, height), top_pixmap, top_pixmap.rect())
        painter.drawPixmap(QRect(width, height * 2, width, height), bottom_pixmap, bottom_pixmap.rect())
        painter.drawPixmap(QRect(0, height, width, height), left_pixmap, left_pixmap.rect())
        painter.drawPixmap(QRect(width * 2, height, width, height), right_pixmap, right_pixmap.rect())
        painter.end()

        Utils.save_pixmap(pixmap)

    @Slot()
    def refresh(self, cells, index):
        self.index = index
        cell = cells[0] if cells else None
        
        self.refresh_probs_tab(cell)
        self.refresh_sensors_tab(cell)
        self.refresh_gene_states_tab(cell)
        self.refresh_interaction_tab(cell, index)

    def refresh_gene_states_tab(self, cell):
        self.clear_chart(self.gene_scores_chart)
        
        if cell is None: #if cell we deselected
            self.gs_table.setRowCount(0) #this deletes all the rows
        else:
            #update the table
            table_data = self.data_tools.get_gs_table_data(cell, self.index) #this is a 2D array
            self.gs_table.setRowCount(len(table_data))

            for row in range(len(table_data)):
                for col in range(len(table_data[row])):
                    item = QTableWidgetItem(table_data[row][col])
                    self.gs_table.setItem(row, col, item)

            #update the scores chart
            scores = self.data_tools.get_gene_scores(self.index)
            num_genes = len(scores)

            x_axis = QtCharts.QBarCategoryAxis()
            categories = list(map(lambda i: str(i), range(num_genes)))
            x_axis.append(categories)
            self.gene_scores_chart.addAxis(x_axis, Qt.AlignBottom)

            y_axis = QtCharts.QValueAxis()
            y_axis.setRange(0.0, 1.0)
            self.gene_scores_chart.addAxis(y_axis, Qt.AlignLeft)

            series = QtCharts.QBarSeries()
            bar_set = QtCharts.QBarSet('Scores')
            for score in scores:
                bar_set.append(score)
                
            series.append(bar_set)
            self.gene_scores_chart.addSeries(series)
            series.attachAxis(y_axis)

    def refresh_probs_tab(self, cell):
        self.probs_chart.removeAllSeries()
        
        if cell is not None:
            info = self.data_tools.get_probs_info_for_cell(cell) #[(sym_label, prob), ...]
            series = QtCharts.QPieSeries()
            for (label, prob) in info:
                qslice = QtCharts.QPieSlice(label, prob)
                series.append(qslice)
                
            series.setLabelsVisible(True)
            self.probs_chart.addSeries(series)

    def refresh_sensors_tab(self, cell):
        sensor_concs = self.data_tools.get_sensor_concs(cell)
        
        pairs = (
            ('Top', self.top_chart),
            ('Right', self.right_chart),
            ('Bottom', self.bottom_chart),
            ('Left', self.left_chart),
        )

        for (loc, chart) in pairs:
            self.clear_chart(chart)

            if cell is not None: #note: cell will be None when an item in the graphics area is de-selected
                num_concs = len(sensor_concs[loc])
                categories = list(map(lambda i: str(i), range(num_concs)))
                x_axis = QtCharts.QBarCategoryAxis()
                x_axis.append(categories)
                chart.addAxis(x_axis, Qt.AlignBottom)

                y_axis = QtCharts.QValueAxis()
                y_axis.setRange(0.0, 1.0)
                chart.addAxis(y_axis, Qt.AlignLeft)

                series = QtCharts.QBarSeries()
                bar_set = QtCharts.QBarSet(loc)
                for sensor_conc in sensor_concs[loc]:
                    bar_set.append(sensor_conc)
                    
                series.append(bar_set)
                chart.addSeries(series)
                series.attachAxis(y_axis)

    def refresh_interaction_tab(self, cell, index):
        deselected = self.interaction_cell and cell == None
        self.interaction_cell = cell
        
        if deselected:
            self.update_interaction_img()

    def update_interaction_img(self):
        if self.interaction_cell is None: #if cell was deselected
            pixmap = QPixmap(0, 0)
            self.interaction_img.setPixmap(pixmap)
            self.interaction_img.resize(0, 0)
        else:
            pixmap_data = self.data_tools.get_interaction_graph(self.index[0], self.index[1], self.interaction_cell)
            pixmap = QPixmap(0, 0)
            pixmap.loadFromData(pixmap_data, format='png')
            size = pixmap.size()
            self.interaction_img.resize(size.width(), size.height())
            self.interaction_img.setPixmap(pixmap)
