from PySide2.QtCore import *
from PySide2.QtGui import *
from PySide2.QtWidgets import *
from PySide2.QtCharts import QtCharts

from Utils import Utils

class CellArea(QWidget):
    def __init__(self, data_tools, *args, **kwargs):
        QWidget.__init__(self, *args, **kwargs)
        self.data_tools = data_tools
        self.cell = None
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
            headers.append('Threshold')

        for i in range(num_bind_sites):
            headers.append('Prod Site {}'.format(i + 1))
            headers.append('Inhib Protein {}'.format(i + 1))
            headers.append('Prod Rate {}'.format(i + 1))

        self.gs_table = QTableWidget(0, len(headers))
        self.gs_table.setHorizontalHeaderLabels(headers)

        #build gene scores chart
        # self.gene_scores_chart = QtCharts.QChart()
        # self.gene_scores_chart.setTitle('Gene Scores')
        # self.gene_scores_view = QtCharts.QChartView(self.gene_scores_chart)
        # self.gene_scores_view.setRenderHint(QPainter.RenderHint.Antialiasing)

        export_button = QPushButton('Export All Steps')
        export_button.setSizePolicy(QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed))
        export_button.clicked.connect(self.export_gs_data)

        #gene scores graph save button
        # save_button = QPushButton('Save')
        # save_button.setSizePolicy(QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed))
        # save_button.clicked.connect(lambda: Utils.save_chart_view(self.gene_scores_view))
        
        #put everything together
        layout.addWidget(self.gs_table)
        layout.addWidget(export_button)
        # layout.addWidget(self.gene_scores_view)
        # layout.addWidget(save_button)
        tab.setLayout(layout)
        
        return tab

    @Slot()
    def export_gs_data(self):
        filename = Utils.run_save_csv_dialog()
        if filename:
            self.data_tools.save_all_gs_table_data(self.cell, self.index, filename)
            QMessageBox.information(self, 'Data Exported', 'Done.')

    def build_interaction_tab(self):
        tab = QWidget()
        layout = QVBoxLayout()
        
        button = QPushButton('Generate')
        button.setSizePolicy(QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed))
        button.clicked.connect(self.update_interaction_img)
        layout.addWidget(button)
        
        scroll_area = QScrollArea()
        inner_layout = QVBoxLayout()
        self.interaction_img = QLabel()
        self.no_interaction_msg = QLabel("No Interaction")
        self.no_interaction_msg.hide()
        inner_layout.addWidget(self.no_interaction_msg)
        inner_layout.addWidget(self.interaction_img)
        scroll_area.setLayout(inner_layout)
        #scroll_area.setWidget(self.interaction_img)
        layout.addWidget(scroll_area)

        save_button = QPushButton('Save')
        save_button.setSizePolicy(QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed))
        save_button.clicked.connect(lambda: Utils.save_pixmap(self.interaction_img.pixmap()))
        save_all_button = QPushButton('Save for all Reg Steps')
        save_all_button.setSizePolicy(QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed))
        save_all_button.clicked.connect(self.save_all_interaction_steps)
        buttons_widget = QWidget()
        buttons_layout = QHBoxLayout()
        buttons_layout.addWidget(save_button)
        buttons_layout.addWidget(save_all_button)
        buttons_widget.setLayout(buttons_layout)
        
        layout.addWidget(buttons_widget)
        tab.setLayout(layout)
        
        return tab

    @Slot()
    def save_all_interaction_steps(self):
        if self.cell:
            path = Utils.run_save_folder_dialog()
            if path:
                self.data_tools.save_all_interaction_graphs(self.index, self.cell, path)
                QMessageBox.information(self, 'Images Saved', 'Done.')

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

        num_locs = 3 + self.data_tools.get_run().max_children
        self.sensor_charts = []
        self.sensor_views = []

        for i in range(num_locs):
            chart = QtCharts.QChart()
            chart.legend().setVisible(False)
            view = QtCharts.QChartView(chart)
            view.setRenderHint(QPainter.RenderHint.Antialiasing)
            chart.resize(*chart_size)

            self.sensor_charts.append(chart)
            self.sensor_views.append(view)

        #these must follow the order given in ProteinPropsMod
        grid_layout.addWidget(self.sensor_views[0], 1, 2) #right
        grid_layout.addWidget(self.sensor_views[1], 0, 1) #top
        grid_layout.addWidget(self.sensor_views[2], 1, 0) #left

        #bottom (children)
        for i in range(3, num_locs):
            grid_layout.addWidget(self.sensor_views[i], 2, i - 3)

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
        num_bottom_locs = self.data_tools.get_run().max_children
        pixmap = QPixmap(max(3, num_bottom_locs) * width, 3 * height)
        pixmap.fill(Qt.white)

        right_pixmap = QPixmap(width, height)
        self.sensor_views[0].render(right_pixmap)
        
        top_pixmap = QPixmap(width, height)
        self.sensor_views[1].render(top_pixmap)

        left_pixmap = QPixmap(width, height)
        self.sensor_views[2].render(left_pixmap)

        bottom_pixmaps = []
        for i in range(num_bottom_locs):
            bpixmap = QPixmap(width, height)
            self.sensor_views[3 + i].render(bpixmap)
            bottom_pixmaps.append(bpixmap)

        painter = QPainter(pixmap)
        painter.drawPixmap(QRect(width * 2, height, width, height), right_pixmap, right_pixmap.rect())
        painter.drawPixmap(QRect(width, 0, width, height), top_pixmap, top_pixmap.rect())
        painter.drawPixmap(QRect(0, height, width, height), left_pixmap, left_pixmap.rect())

        for i in range(len(bottom_pixmaps)):
            painter.drawPixmap(QRect(i * width, height * 2, width, height), bottom_pixmaps[i], bottom_pixmaps[i].rect())
            
        painter.end()
        Utils.save_pixmap(pixmap)

    @Slot()
    def refresh(self, cells, index):
        self.index = index
        cell = cells[0] if cells else None
        deselected = self.cell and cell == None
        self.cell = cell
        
        self.refresh_probs_tab()
        self.refresh_sensors_tab()
        self.refresh_gene_states_tab()
        self.refresh_interaction_tab(deselected)

    def refresh_gene_states_tab(self):
        #self.clear_chart(self.gene_scores_chart)
        
        if self.cell is None: #if self.cell we deselected
            self.gs_table.setRowCount(0) #this deletes all the rows
        else:
            #update the table
            table_data = self.data_tools.get_gs_table_data(self.cell, self.index) #this is a 2D array
            self.gs_table.setRowCount(len(table_data))

            for row in range(len(table_data)):
                for col in range(len(table_data[row])):
                    item = QTableWidgetItem(table_data[row][col])
                    self.gs_table.setItem(row, col, item)

            #update the scores chart
            #scores = self.data_tools.get_gene_scores(self.index)
            #num_genes = len(scores)

            # x_axis = QtCharts.QBarCategoryAxis()
            # categories = list(map(lambda i: str(i), range(num_genes)))
            # x_axis.append(categories)
            # self.gene_scores_chart.addAxis(x_axis, Qt.AlignBottom)

            # y_axis = QtCharts.QValueAxis()
            # #y_axis.setRange(0.0, 1.0)
            # self.gene_scores_chart.addAxis(y_axis, Qt.AlignLeft)

            # series = QtCharts.QBarSeries()
            # bar_set = QtCharts.QBarSet('Scores')
            # for score in scores:
            #     bar_set.append(score)
                
            # series.append(bar_set)
            # self.gene_scores_chart.addSeries(series)
            # series.attachAxis(y_axis)

    def refresh_probs_tab(self):
        self.probs_chart.removeAllSeries()
        
        if self.cell is not None:
            info = self.data_tools.get_probs_info_for_cell(self.cell) #[(sym_label, prob), ...]
            series = QtCharts.QPieSeries()
            for (label, prob) in info:
                qslice = QtCharts.QPieSlice(label, prob)
                series.append(qslice)
                
            series.setLabelsVisible(True)
            self.probs_chart.addSeries(series)

    def refresh_sensors_tab(self):
        sensor_concs = self.data_tools.get_sensor_concs(self.cell)

        for loc in range(len(sensor_concs)):
            chart = self.sensor_charts[loc]
            self.clear_chart(chart)

            if self.cell is not None: #note: self.cell will be None when an item in the graphics area is de-selected
                num_concs = len(sensor_concs[loc])
                categories = list(map(lambda i: str(i), range(num_concs)))
                x_axis = QtCharts.QBarCategoryAxis()
                x_axis.append(categories)
                chart.addAxis(x_axis, Qt.AlignBottom)

                y_axis = QtCharts.QValueAxis()
                y_axis.setRange(0.0, 1.0)
                chart.addAxis(y_axis, Qt.AlignLeft)

                series = QtCharts.QBarSeries()
                bar_set = QtCharts.QBarSet(str(loc))
                for sensor_conc in sensor_concs[loc]:
                    bar_set.append(sensor_conc)
                    
                series.append(bar_set)
                chart.addSeries(series)
                series.attachAxis(y_axis)

    def refresh_interaction_tab(self, deselected):
        if deselected:
            self.update_interaction_img()

    @Slot()
    def update_interaction_img(self):
        if self.cell is None: #if cell was deselected
            pixmap = QPixmap(0, 0)
            self.interaction_img.setPixmap(pixmap)
            self.interaction_img.resize(0, 0)
            self.interaction_img.hide()
            self.no_interaction_msg.hide()
        else:
            pixmap_data = self.data_tools.get_interaction_graph(self.index, self.cell)
            if pixmap_data is None:
                self.no_interaction_msg.show()
                self.interaction_img.hide()
            else:
                self.no_interaction_msg.hide()
                pixmap = QPixmap(0, 0)
                pixmap.loadFromData(pixmap_data, format='png')
                size = pixmap.size()
                self.interaction_img.resize(size.width(), size.height())
                self.interaction_img.setPixmap(pixmap)
                self.interaction_img.show()
