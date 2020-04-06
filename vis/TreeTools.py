from graphviz import Digraph
from DataTools import DataTools

class TreeTools():
    def __init__(self, data):
        self.data = data
    
    def gen_image(self, data, index):
        tree = data.get_tree(index)
        graph = Digraph(format='png')
        self.build_graph(graph, None, tree.root)

        return graph.pipe()

    def build_graph(self, graph, parent_name, cell):
        name = str(id(cell)) #use pointer value as unique identifier
        label = str(cell.sym.val) if cell.sym != None else '_'
        graph.node(name, label)

        #add edges from parent (if one exists)
        if parent_name is not None:
            graph.edge(parent_name, name)

        for child in cell.children:
            self.build_graph(graph, name, child)
    
