from PySide2.QtCore import *
from PySide2.QtGui import *
from PySide2.QtWidgets import *
from PySide2.QtCharts import *

from Utils import Utils

class FitnessArea(QWidget):
    def __init__(self, data_tools, *args, **kwargs):
        QWidget.__init__(self, *args, **kwargs)
        self.data_tools = data_tools
        self.run = self.data_tools.get_run()

        self.chart_view = QtCharts.QChartView()
        self.gen_button = QPushButton('Generate')
        self.gen_button.setSizePolicy(QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed))
        self.gen_button.clicked.connect(self.update_chart)
        
        save_button = QPushButton('Save')
        save_button.setSizePolicy(QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed))
        save_button.clicked.connect(lambda: Utils.save_chart_view(self.chart_view))

        self.show_best_checkbox = QCheckBox("Show Best")
        self.show_best_checkbox.setCheckState(Qt.Checked)
        self.show_gbest_checkbox = QCheckBox("Show Gen Best")
        self.show_gbest_checkbox.setCheckState(Qt.Checked)
        self.show_gbest_avg_checkbox = QCheckBox("Show Gen Avg")
        self.show_gbest_avg_checkbox.setCheckState(Qt.Checked)

        self.show_best_checkbox
        
        layout = QVBoxLayout()
        checks_group_widget = QWidget()
        checks_vbox = QVBoxLayout()
        checks_vbox.addWidget(self.show_best_checkbox)
        checks_vbox.addWidget(self.show_gbest_checkbox)
        checks_vbox.addWidget(self.show_gbest_avg_checkbox)
        checks_vbox.addWidget(self.gen_button)
        checks_group_widget.setLayout(checks_vbox)
        
        layout.addWidget(checks_group_widget)
        layout.addWidget(self.chart_view)
        layout.addWidget(save_button)

        self.setLayout(layout)

    @Slot()
    def update_chart(self):
        bests, gen_bests, gen_avgs = self.data_tools.get_best_fitnesses()
        chart = QtCharts.QChart()

        x_axis = QtCharts.QValueAxis()
        x_axis.setRange(0, self.run.ea_steps)
        #x_axis.setTickInterval(self.run.step_range.step)
        x_axis.setTickCount(self.run.ea_steps // self.run.step_range.step + 1)        
        x_axis.setLabelFormat('%0.0f') #display as integer
        x_axis.setTitleText('EA Step')

        y_axis = QtCharts.QValueAxis()
        y_axis.setRange(0.0, 1.0)
        y_axis.setTitleText('Fitness')

        chart.addAxis(x_axis, Qt.AlignBottom)
        chart.addAxis(y_axis, Qt.AlignLeft)

        if self.show_best_checkbox.isChecked():
            best_series = QtCharts.QLineSeries()
            best_series.setName('Best Fitness')
            chart.addSeries(best_series)
            for ea_step in range(0, self.run.ea_steps + 1):
                best_series.append(QPointF(ea_step, bests[ea_step]))
                
            best_series.attachAxis(x_axis)
            best_series.attachAxis(y_axis)

        if self.show_gbest_checkbox.isChecked():
            gen_best_series = QtCharts.QLineSeries()
            gen_best_series.setName('Gen Best Fitness')
            chart.addSeries(gen_best_series)
            for ea_step in range(0, self.run.ea_steps + 1):
                gen_best_series.append(QPointF(ea_step, gen_bests[ea_step]))
                
            gen_best_series.attachAxis(x_axis)
            gen_best_series.attachAxis(y_axis)

        if self.show_gbest_avg_checkbox.isChecked():
            avg_series = QtCharts.QLineSeries()
            avg_series.setName('Gen Avg')
            chart.addSeries(avg_series)
            for ea_step in range(0, self.run.ea_steps + 1):
                avg_series.append(QPointF(ea_step, gen_avgs[ea_step]))
                
            avg_series.attachAxis(x_axis)
            avg_series.attachAxis(y_axis)

        chart.legend().setVisible(True)
        chart.legend().setAlignment(Qt.AlignBottom)

        self.chart_view.setChart(chart)
        self.chart_view.setRenderHint(QPainter.RenderHint.Antialiasing)
