from PySide2.QtCore import *
from PySide2.QtGui import *
from PySide2.QtWidgets import *
from PySide2.QtCharts import *

from CustomStyle import CustomStyle

class FitnessBreakdownArea(QWidget):
    def __init__(self, data_tools, *args, **kwargs):
        QWidget.__init__(self, *args, **kwargs)
        self.data_tools = data_tools
        self.run = self.data_tools.get_run()
        rows = self.data_tools.get_pop_fitness_breakdown()

        layout = QVBoxLayout()
        tabs = QTabWidget(self)
        table_tab = self.build_table_tab(rows)
        graph_tab = self.build_graph_tab(rows)

        tabs.addTab(table_tab, "Table")
        tabs.addTab(graph_tab, "Graph")

        layout.addWidget(tabs)
        self.setLayout(layout)

    def build_graph_tab(self, rows):
        tab = QWidget()
        layout = QVBoxLayout()

        graph_series = {}
        checks = []
        for i in range(len(rows)):
            if i == 0: #header row
                for label in rows[i][1:]: #skip EA Step column
                    series = QtCharts.QLineSeries()
                    series.setName(label)
                    graph_series[label] = series

                    cbox = QCheckBox(label)
                    checks.append(cbox)

            else:
                ea_step = rows[i][0]
                for j in range(1, len(rows[i])): #skip EA Step column
                    label = rows[0][j]
                    val = rows[i][j]
                    series = graph_series[label]
                    series.append(QPointF(float(ea_step), float(val)))


        for cbox in checks:
            layout.addWidget(cbox)
            cbox.stateChanged.connect(lambda new_val: self.update_graph(checks, graph_series))
            
        chart_view = self.build_graph(checks, graph_series)

        layout.addWidget(chart_view)
        tab.setLayout(layout)

        return tab

    @Slot()
    def update_graph(self, checks, graph_series):
        for cbox in checks:
            label = cbox.text()
            series = graph_series[label]
            series.setVisible(cbox.isChecked())

    def build_graph(self, checks, graph_series):
        chart = QtCharts.QChart()

        x_axis = QtCharts.QValueAxis()
        x_axis.setRange(0, self.run.ea_steps)
        x_axis.setTickCount(self.run.ea_steps // self.run.step_range.step + 1)        
        x_axis.setLabelFormat('%0.0f') #display as integer
        x_axis.setTitleText('EA Step')

        y_axis = QtCharts.QValueAxis()
        y_axis.setRange(0.0, 1.0)
        y_axis.setTitleText('Fitness')

        chart.addAxis(x_axis, Qt.AlignBottom)
        chart.addAxis(y_axis, Qt.AlignLeft)

        for cbox in checks:
            label = cbox.text()
            series = graph_series[label]
            series.setVisible(cbox.isChecked())
            chart.addSeries(series)
            series.attachAxis(x_axis)
            series.attachAxis(y_axis)

        chart.legend().setVisible(True)
        chart.legend().setAlignment(Qt.AlignBottom)

        chart_view = QtCharts.QChartView()
        chart_view.setChart(chart)
        chart_view.setRenderHint(QPainter.RenderHint.Antialiasing)

        return chart_view
        
    def build_table_tab(self, rows):
        tab = QWidget()
        layout = QVBoxLayout()
        label = QLabel("Averages for each ea_step:")
        layout.addWidget(label)
        table = self.build_table(rows)
        layout.addWidget(table)
        tab.setLayout(layout)

        return tab
        
    def build_table(self, rows):
        headers = rows[0]
        table = QTableWidget(len(rows) - 1, len(headers))
        table.setHorizontalHeaderLabels(headers)
        table.verticalHeader().setVisible(False)

        widths = list(map(len, rows[0])) #start with widths of headers
        num_cols = len(rows[0]) #rows[0] is headers
        
        for i in range(1, len(rows)): #start at 1 to skip over header row
            for j in range(num_cols):
                if j < len(rows[i]): #indiv genomes may grow over time...
                    item = QTableWidgetItem(rows[i][j])
                    table.setItem(i - 1, j, item)
                else:
                    item = QTableWidgetItem('')
                    table.setItem(i - 1, j, item)

                widths[j] = max(widths[j], len(item.text()))

        for col in range(num_cols):
            table.setColumnWidth(col, widths[col] * 7)

        return table
