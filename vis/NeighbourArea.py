from PySide2.QtCore import *
from PySide2.QtGui import *
from PySide2.QtWidgets import *
from PySide2.QtCharts import *

from Utils import Utils

class NeighbourArea(QWidget):
    def __init__(self, data_tools, initial_index, *args, **kwargs):
        QWidget.__init__(self, *args, **kwargs)
        self.data_tools = data_tools
        self.index = initial_index

        scroll_area = QScrollArea()
        inner_layout = QVBoxLayout()
        self.img = QLabel()
        inner_layout.addWidget(self.img)
        scroll_area.setLayout(inner_layout)
        
        self.gen_button = QPushButton('Generate')
        self.gen_button.setSizePolicy(QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed))
        self.gen_button.clicked.connect(self.update_img)
        
        save_button = QPushButton('Save')
        save_button.setSizePolicy(QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed))
        save_button.clicked.connect(lambda: Utils.save_pixmap(self.img.pixmap()))

        layout = QVBoxLayout()
        layout.addWidget(self.gen_button)
        layout.addWidget(self.img)
        layout.addWidget(save_button)
        self.setLayout(layout)

    @Slot()
    def update(self, index):
        self.index = index
        self.img.hide()

    @Slot()
    def update_img(self):
        pixmap_data = self.data_tools.get_neighbour_graph(self.index)
        if pixmap_data is not None:
            pixmap = QPixmap(0, 0)
            pixmap.loadFromData(pixmap_data, format='png')
            size = pixmap.size()
            self.img.resize(size.width(), size.height())
            self.img.setPixmap(pixmap)
            self.img.show()
