from PySide2.QtCore import *
from PySide2.QtGui import *
from PySide2.QtWidgets import *
from PySide2.QtCharts import *

class FitnessArea(QWidget):
    def __init__(self, data_tools, *args, **kwargs):
        QWidget.__init__(self, *args, **kwargs)
        self.data_tools = data_tools
        self.run = self.data_tools.get_run()

        self.chart_view = QtCharts.QChartView()
        self.gen_button = QPushButton('Generate')
        self.gen_button.setSizePolicy(QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed))
        self.gen_button.clicked.connect(self.update_chart)

        layout = QVBoxLayout()
        layout.addWidget(self.gen_button)
        layout.addWidget(self.chart_view)

        self.setLayout(layout)

    @Slot()
    def update_chart(self):
        bests, gen_bests = self.data_tools.get_best_fitnesses()
        chart = QtCharts.QChart()

        x_axis = QtCharts.QValueAxis()
        x_axis.setRange(0, self.run.ea_steps)
        x_axis.setTickInterval(self.run.step_range.step)
        x_axis.setTickCount(self.run.ea_steps // self.run.step_range.step + 1)        
        x_axis.setLabelFormat('%0.0f') #display as integer
        x_axis.setTitleText('EA Step')

        y_axis = QtCharts.QValueAxis()
        y_axis.setRange(0.0, 2.0)
        y_axis.setTitleText('Fitness')

        chart.addAxis(x_axis, Qt.AlignBottom)
        chart.addAxis(y_axis, Qt.AlignLeft)

        best_series = QtCharts.QLineSeries()
        best_series.setName('Best Fitness')
        gen_best_series = QtCharts.QLineSeries()
        gen_best_series.setName('Gen Best Fitness')
        chart.addSeries(gen_best_series)
        chart.addSeries(best_series)

        i = 0
        for ea_step in range(0, self.run.ea_steps + 1, self.run.step_range.step):
            best_series.append(QPointF(ea_step, bests[i]))
            gen_best_series.append(QPointF(ea_step, gen_bests[i]))
            i += 1

        best_series.attachAxis(x_axis)
        best_series.attachAxis(y_axis)
        gen_best_series.attachAxis(x_axis)
        gen_best_series.attachAxis(y_axis)

        chart.legend().setVisible(True)
        chart.legend().setAlignment(Qt.AlignBottom)

        self.chart_view.setChart(chart)
        self.chart_view.setRenderHint(QPainter.RenderHint.Antialiasing)
        
        self.gen_button.hide()        
