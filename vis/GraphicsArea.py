from PySide2.QtCore import *
from PySide2.QtGui import *
from PySide2.QtWidgets import *

from Utils import Utils
from CustomGraphicsPixmapItem import CustomGraphicsPixmapItem

class GraphicsArea(QWidget):
    selectionChanged = Signal(list)
    
    def __init__(self, data_tools, tree_tools, initial_index, *args, **kwargs):
        QWidget.__init__(self, *args, **kwargs)
        self.data_tools = data_tools
        layout = QVBoxLayout()

        self.view = CustomGraphicsView(tree_tools, initial_index)
        save_button = QPushButton('Save')
        save_button.setSizePolicy(QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed))
        self.fitness_label = QLabel('Fitness: ')
        self.fitness_label.setTextInteractionFlags(Qt.TextSelectableByMouse)
        self.code_label = QLabel('Code: ')
        self.code_label.setTextInteractionFlags(Qt.TextSelectableByMouse)
        
        save_button.clicked.connect(lambda: Utils.save_graphics_view(self.view))
        self.view.selectionChanged.connect(self.handle_selectionChanged)

        layout.addWidget(self.view)
        layout.addWidget(self.fitness_label)
        layout.addWidget(self.code_label)
        layout.addWidget(save_button)

        self.setLayout(layout)
        self.update_info_labels(initial_index)

    @Slot()
    def refresh(self, index, checked_info=[]):
        self.view.refresh(index, checked_info)
        self.update_info_labels(index)

    @Slot()
    def handle_selectionChanged(self, cells):
        self.selectionChanged.emit(cells)

    def disconnect_signals(self):
        self.view.disconnect_signals()

    def update_info_labels(self, index):
        fitness = self.data_tools.get_indiv_fitness(index)
        code = self.data_tools.get_indiv_code(index)
        self.fitness_label.setText('Fitness: {:.2f}'.format(fitness))
        self.code_label.setText('Code: {}'.format(code))

class CustomGraphicsView(QGraphicsView):
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
        sel_items = self.scene.selectedItems()
        sel_cell = sel_items[0].cell if sel_items else None
        
        self.tree_tools.clear_scene(self.scene) #clear any previous tree
        self.tree_tools.draw_scene(self.scene, index, checked_info)

        if sel_cell is not None:
            self.reselect_cell(sel_cell)

    @Slot()
    def handle_selectionChanged(self):
        cells = list(map(lambda item: item.cell, self.scene.selectedItems()))
        self.selectionChanged.emit(cells)

    def disconnect_signals(self):
        self.scene.selectionChanged.disconnect(self.handle_selectionChanged)

    def reselect_cell(self, cell):
        sel_id = cell.id
        items = self.scene.items()

        i = 0
        found = False
        while not found and i < len(items):
            if isinstance(items[i], CustomGraphicsPixmapItem) and items[i].cell.id == sel_id:
                items[i].setSelected(True)
                found = True

            i += 1
