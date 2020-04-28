from PySide2.QtGui import *
from PySide2.QtWidgets import *
from PySide2.QtCore import *
from PySide2.QtCharts import *

from DataTools import DataTools
from DrawTree import DrawTree
from TreeLayout import TreeLayout
from CustomGraphicsPixmapItem import CustomGraphicsPixmapItem
from SettingsArea import Settings

class TreeTools():
    node_width = 200
    node_height = 200
    node_space = 20

    def __init__(self, data_tools):
        self.data_tools = data_tools
        self.run = data_tools.get_run()

    def clear_scene(self, scene):
        for item in scene.items():
            scene.removeItem(item)

    #call this one!
    def draw_scene(self, scene, index, checked_info):
        cell_tree = self.data_tools.get_tree(index)
        root_cell = self.data_tools.get_root_cell(cell_tree)
        
        tree = DrawTree(root_cell, self.data_tools)
        tree = TreeLayout.buchheim(tree)

        self.draw_tree(scene, tree, 0, checked_info)
        self.draw_edges(scene, tree, 0)
                  
    def draw_tree(self, scene, node, depth, checked_info):
        pixmap = self.build_conc_graph(node.cell, checked_info)
        item = CustomGraphicsPixmapItem(pixmap, node.cell)
        item.setPos(QPointF(
            node.x * TreeTools.node_width + node.x * TreeTools.node_space,
            depth * TreeTools.node_height + node.y * TreeTools.node_space
        ))
        item.setFlags(QGraphicsItem.ItemIsSelectable)
        scene.addItem(item)
        
        for child in node.children:
            self.draw_tree(scene, child, depth + 1, checked_info)

    def draw_edges(self, scene, node, depth):
        pen = QPen(Qt.black, 1)
        for child in node.children:
            line = QGraphicsLineItem(
                node.x * TreeTools.node_width + (TreeTools.node_width / 2) + node.x * TreeTools.node_space,
                (depth + 1) * TreeTools.node_height + (depth) * TreeTools.node_space - 5,

                child.x * TreeTools.node_width + (TreeTools.node_width / 2) + child.x * TreeTools.node_space,
                (depth + 1) * TreeTools.node_height + (depth + 1) * TreeTools.node_space + 5
            )
            line.setPen(pen)
            scene.addItem(line)
            
            self.draw_edges(scene, child, depth + 1)
    
    def build_conc_graph(self, cell, checked_info):
        chart = QtCharts.QChart() #note: we're not going to show this widget in the GUI, it's just for generating charts
        num_concs = self.data_tools.get_num_genes(cell)

        #build bar_series
        bar_series = QtCharts.QBarSeries()
        for props, colour in checked_info:
            protein = self.data_tools.get_protein(cell, props)
            if protein is not None:
                concs = protein.concs
            else:
                concs = []

            props_str = self.data_tools.get_props_str(props)
            bar_set = QtCharts.QBarSet(props_str)
            bar_set.setColor(colour)
            for conc in concs:
                bar_set.append(conc)
            bar_series.append(bar_set)
            
        #line_series
        threshold_series = self.get_threshold_series(num_concs)

        #add series to chart
        chart.addSeries(bar_series)
        for line_series in threshold_series:
            chart.addSeries(line_series)

        #create x axis and add to chart
        categories = list(map(str, range(num_concs)))
        x_axis = QtCharts.QBarCategoryAxis()
        x_axis.append(categories)
        chart.addAxis(x_axis, Qt.AlignBottom)

        #attach x axis to series
        bar_series.attachAxis(x_axis)
        for line_series in threshold_series:
            line_series.attachAxis(x_axis)

        #create y axis and add to chart
        y_axis = QtCharts.QValueAxis()
        chart.addAxis(y_axis, Qt.AlignLeft)

        #attach y axis to series
        bar_series.attachAxis(y_axis)
        for line_series in threshold_series:
            line_series.attachAxis(y_axis)
            
        y_axis.setRange(0.0, 1.0)

        #title and legend
        cell_label = str(cell.sym.val) if cell.sym != None else '(None)'
        chart.setTitle(cell_label)
        chart.legend().setVisible(False)

        #view
        chartView = QtCharts.QChartView(chart)
        chartView.setRenderHint(QPainter.RenderHint.Antialiasing)
        chartView.resize(TreeTools.node_width, TreeTools.node_height)
        pixmap = QPixmap(chartView.size())
        chartView.render(pixmap)

        return pixmap

    def get_threshold_series(self, num_concs):
        series = []
        pairs = (
            (Settings.show_protein_deletion_threshold, self.run.protein_deletion_threshold),
            (Settings.show_cell_division_threshold, self.run.cell_division_threshold),
            (Settings.show_sensor_reinforcement_threshold, self.run.sensor_reinforcement_threshold),
            (Settings.show_sym_prob_threshold, self.run.sym_prob_threshold),
        )
        for (setting, threshold) in pairs:
            if setting:
                cur_series = self.get_single_threshold_series(threshold, num_concs)
                series.append(cur_series)

        return series

    def get_single_threshold_series(self, threshold, num_concs):
        series = QtCharts.QLineSeries()
        pen = QPen(Qt.DashLine)
        series.setPen(pen)
        #just make a straight line - overshoot the boundaries so the line goes all the way across the graph
        for i in range(-1, num_concs + 1):
            series.append(QPointF(i, threshold))

        return series
            
