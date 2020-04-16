from PySide2.QtWidgets import QGraphicsPixmapItem

class CustomGraphicsPixmapItem(QGraphicsPixmapItem):
    def __init__(self, pixmap, cell):
        QGraphicsPixmapItem.__init__(self, pixmap)
        self.cell = cell
