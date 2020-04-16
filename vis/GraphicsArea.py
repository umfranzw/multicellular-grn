from PySide2.QtCore import *
from PySide2.QtGui import *
from PySide2.QtWidgets import *

class GraphicsArea(QGraphicsView):
    selectionChanged = Signal(list)
    
    def __init__(self, tree_tools, initial_index, *args, **kwargs):
        QGraphicsView.__init__(self, *args, **kwargs)
        self.scene = QGraphicsScene()
        #self.scene.setBackgroundBrush(QBrush(Qt.black))
        self.setScene(self.scene)
        
        self.tree_tools = tree_tools

        self.scene.selectionChanged.connect(self.handle_selectionChanged)

        self.refresh(initial_index)

    @Slot()
    def refresh(self, index, checked_info=[]):
        self.tree_tools.clear_scene(self.scene) #clear any previous tree
        self.tree_tools.draw_scene(self.scene, index, checked_info)   

    @Slot()
    def handle_selectionChanged(self):
        cells = list(map(lambda item: item.cell, self.scene.selectedItems()))
        self.selectionChanged.emit(cells)
