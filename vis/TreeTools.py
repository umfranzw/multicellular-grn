from PySide2.QtGui import QImage, QPainter, QPixmap, QPen
from PySide2.QtWidgets import QGraphicsItem, QGraphicsLineItem
from DataTools import DataTools
from PySide2.QtCharts import QtCharts
from PySide2.QtCore import Qt, QPointF
from DrawTree import DrawTree
from TreeLayout import TreeLayout
from CustomGraphicsPixmapItem import CustomGraphicsPixmapItem

class TreeTools():
    node_width = 200
    node_height = 200
    node_space = 50

    def __init__(self, data_tools):
        self.data_tools = data_tools

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
        item.setOffset(QPointF(node.x * TreeTools.node_width, depth * TreeTools.node_height))
        item.setFlags(QGraphicsItem.ItemIsSelectable)
        scene.addItem(item)
        
        for child in node.children:
            draw_tree(scene, child, depth + 1, checked_info)

    def draw_edges(self, scene, node, depth):
        pen = QPen(Qt.black, 1)
        for child in node.children:
            line = QGraphicsLineItem(node.x * TreeTools.node_width + (TreeTools.node_space / 2), depth * TreeTools.node_height + (TreeTools.node_space / 2),
                 child.x * TreeTools.node_width + (TreeTools.node_space / 2), (depth + 1) * TreeTools.node_height + (TreeTools.node_space / 2))
            line.setPen(pen)
            scene.addItem(line)
            
            self.draw_edges(scene, child, depth + 1)
    
    def build_conc_graph(self, cell, checked_info):
        chart = QtCharts.QChart() #note: we're not going to show this widget in the GUI, it's just for generating charts
        cell_label = str(cell.sym.val) if cell.sym != None else '_'
        chart.setTitle(cell_label)

        num_concs = self.data_tools.get_num_genes(cell)
        categories = list(map(lambda i: str(i), range(num_concs)))
        x_axis = QtCharts.QBarCategoryAxis()
        x_axis.append(categories)
        chart.addAxis(x_axis, Qt.AlignBottom)

        y_axis = QtCharts.QValueAxis()
        y_axis.setRange(0.0, 1.0)
        chart.addAxis(y_axis, Qt.AlignLeft)
        
        series = QtCharts.QBarSeries()
        for props, colour in checked_info:
            protein = self.data_tools.get_protein(cell, props)
            if protein is not None:
                concs = protein.concs
            else:
                concs = []

            props_str = self.data_tools.get_props_str(props)
            bar_set = QtCharts.QBarSet(props_str)
            bar_set.setColor(colour)
            bar_set.append(concs)
            series.append(bar_set)

        
        chart.addSeries(series)
        series.attachAxis(y_axis)
        chart.legend().setVisible(False)
        
        chartView = QtCharts.QChartView(chart)
        chartView.setRenderHint(QPainter.RenderHint.Antialiasing)
        chartView.resize(200, 200)
        pixmap = QPixmap(chartView.size())
        chartView.render(pixmap)

        return pixmap
