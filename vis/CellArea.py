from PySide2.QtCore import *
from PySide2.QtGui import *
from PySide2.QtWidgets import *
from PySide2.QtCharts import QtCharts

class CellArea(QWidget):
    def __init__(self, data_tools, *args, **kwargs):
        QWidget.__init__(self, *args, **kwargs)
        self.data_tools = data_tools
        
        layout = QVBoxLayout()
        tabs = QTabWidget(self)
        probs_tab = self.build_probs_tab()
        sensors_tab = self.build_sensors_tab()
        gene_states_tab = QWidget()

        tabs.addTab(probs_tab, "Sym Probs")
        tabs.addTab(sensors_tab, "Sensors")
        tabs.addTab(gene_states_tab, "Gene States")

        layout.addWidget(tabs)

        self.setLayout(layout)

    def build_probs_tab(self):
        tab = QWidget()
        layout = QVBoxLayout()
        self.probs_chart = QtCharts.QChart()
        self.probs_chart.legend().setVisible(False)
        self.probs_chart.setTitle('Symbol Probs')

        chart_view = QtCharts.QChartView(self.probs_chart)
        chart_view.setRenderHint(QPainter.RenderHint.Antialiasing)
        self.probs_chart.resize(300, 300)
        
        layout.addWidget(chart_view)
        tab.setLayout(layout)

        scroll_area = QScrollArea()
        scroll_area.setWidget(tab)
        
        return scroll_area

    def clear_sensors_chart(self, chart):
        chart.removeAllSeries()
        for axis in chart.axes(Qt.Horizontal | Qt.Vertical):
            chart.removeAxis(axis)

    def refresh_sensors_tab(self, cell):
        sensor_concs = self.data_tools.get_sensor_concs(cell)
        
        pairs = (
            ('Top', self.top_chart),
            ('Right', self.right_chart),
            ('Bottom', self.bottom_chart),
            ('Left', self.left_chart),
        )

        for (loc, chart) in pairs:
            self.clear_sensors_chart(chart)

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
                bar_set.append(sensor_concs[loc])
                series.append(bar_set)
                chart.addSeries(series)
                series.attachAxis(y_axis)

    def build_sensors_tab(self):
        tab = QWidget()
        layout = QGridLayout()
        chart_size = (200, 200)

        self.top_chart = QtCharts.QChart()
        self.top_chart.legend().setVisible(False)
        top_view = QtCharts.QChartView(self.top_chart)
        top_view.setRenderHint(QPainter.RenderHint.Antialiasing)
        self.top_chart.resize(*chart_size)

        self.right_chart = QtCharts.QChart()
        self.right_chart.legend().setVisible(False)
        right_view = QtCharts.QChartView(self.right_chart)
        right_view.setRenderHint(QPainter.RenderHint.Antialiasing)
        self.right_chart.resize(*chart_size)

        self.bottom_chart = QtCharts.QChart()
        self.bottom_chart.legend().setVisible(False)
        bottom_view = QtCharts.QChartView(self.bottom_chart)
        bottom_view.setRenderHint(QPainter.RenderHint.Antialiasing)
        self.bottom_chart.resize(*chart_size)

        self.left_chart = QtCharts.QChart()
        self.left_chart.legend().setVisible(False)
        left_view = QtCharts.QChartView(self.left_chart)
        left_view.setRenderHint(QPainter.RenderHint.Antialiasing)
        self.left_chart.resize(*chart_size)

        layout.addWidget(top_view, 0, 1)
        layout.addWidget(right_view, 1, 2)
        layout.addWidget(bottom_view, 2, 1)
        layout.addWidget(left_view, 1, 0)
        
        tab.setLayout(layout)
        scroll_area = QScrollArea()
        scroll_area.setWidget(tab)

        return scroll_area

    @Slot()
    def refresh(self, cells):
        cell = cells[0] if cells else None
        self.refresh_probs_tab(cell)
        self.refresh_sensors_tab(cell)

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
            
