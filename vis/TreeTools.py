from graphviz import Graph
from PySide2.QtGui import QImage
from DataTools import DataTools
from PySide2.QtCharts import QtCharts
from PySide2.QtCore import Qt
from PySide2.QtGui import QPainter, QPixmap
import tempfile
import shutil

class TreeTools():
    def __init__(self, data):
        self.data = data
        self.img_path = tempfile.mkdtemp()

    def close(self):
        shutil.rmtree(self.img_path)
    
    def gen_image(self, data, index, props_list):
        tree = data.get_tree(index)
        graph = Graph(format='png')
        self.build_graph(graph, None, tree.root, props_list)

        return QImage.fromData(graph.pipe(), format='png')

    #props is the props of the protein to graph (if it's present in cell)
    def build_graph(self, graph, parent_name, cell, checked_info):
        name = str(cell.jl_value) #use Julia pointer value as unique identifier
        label = str(cell.sym.val) if cell.sym != None else '_'

        pixmap = self.build_conc_graph(cell, checked_info)
        filename = '{}/{}.png'.format(self.img_path, name)
        pixmap.save(filename, format='png')
        graph.node(name, label, image=filename, shape='box', labelloc='b')

        #add edges from parent (if one exists)
        if parent_name is not None:
            graph.edge(parent_name, name)

        for child in cell.children:
            self.build_graph(graph, name, child)

    def build_conc_graph(self, cell, checked_info):
        chart = QtCharts.QChart() #note: we're not going to show this widget in the GUI, it's just for generating charts
        series = QtCharts.QBarSeries()
        for props, colour in checked_info:
            protein = self.data.get_protein(cell, props)
            if protein is not None:
                concs = protein.concs
            else:
                concs = []

            props_str = self.data.get_props_str(props)
            bar_set = QtCharts.QBarSet(props_str)
            bar_set.setColor(colour)
            bar_set.append(concs)
            series.append(bar_set)

        chart.addSeries(series)
        num_concs = series.barSets()[0].count() if series.barSets() else 0
        categories = list(map(lambda i: str(i), range(num_concs)))
        x_axis = QtCharts.QBarCategoryAxis()
        x_axis.append(categories)
        chart.addAxis(x_axis, Qt.AlignBottom)

        y_axis = QtCharts.QValueAxis()
        y_axis.setRange(0, 1)
        chart.addAxis(y_axis, Qt.AlignLeft)
        
        chart.legend().setVisible(False)
        
        chartView = QtCharts.QChartView(chart)
        chartView.setRenderHint(QPainter.RenderHint.Antialiasing)
        chartView.resize(200, 200)
        pixmap = QPixmap(chartView.size())
        chartView.render(pixmap)

        return pixmap
