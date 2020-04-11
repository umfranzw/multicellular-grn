from graphviz import Graph
from PySide2.QtGui import QImage
from DataTools import DataTools
import pyqtgraph as pg
import pyqtgraph.exporters
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
    def build_graph(self, graph, parent_name, cell, props_list):
        name = str(cell.jl_value) #use Julia pointer value as unique identifier
        label = str(cell.sym.val) if cell.sym != None else '_'

        qimg = self.build_conc_graph(cell, props_list)
        filename = '{}/{}.png'.format(self.img_path, name)
        qimg.save(filename, format='png')
        graph.node(name, label, image=filename, shape='box')

        #add edges from parent (if one exists)
        if parent_name is not None:
            graph.edge(parent_name, name)

        for child in cell.children:
            self.build_graph(graph, name, child)

    def build_conc_graph(self, cell, props_list):
        plot = pg.PlotWidget() #note: we're not going to show this widget in the GUI, it's just for generating plots
        max_xs_len = 0

        for props in props_list:
            protein = self.data.get_protein(cell, props)
            if protein is not None:
                concs = protein.concs
            else:
                concs = []
            xs = list(range(len(concs)))
            max_xs_len = max(max_xs_len, len(xs))

            bars = pg.BarGraphItem(x=xs, height=concs, width=0.3, brush=pyqtgraph.intColor(hash(props)))
            plot.addItem(bars)
            
        plot.setXRange(0, max_xs_len)
        plot.setYRange(0, 1)

        exporter = pg.exporters.ImageExporter(plot.plotItem)
        exporter.parameters()['width'] = 100 #heigh will be auto-adjusted to match
        qimg = exporter.export(toBytes=True) #generates a QImage

        return qimg

        

        
