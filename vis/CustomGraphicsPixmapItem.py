from PySide2.QtWidgets import *
from PySide2.QtGui import *
from PySide2.QtCore import *

class CustomGraphicsPixmapItem(QGraphicsPixmapItem):
    def __init__(self, pixmap, cell, *args, **kwargs):
        painter = QPainter(pixmap)
        painter.setPen(Qt.black)
        rect = QRectF(0, pixmap.height() - 30, pixmap.width(), 20)
        painter.drawText(rect, 'Age: {}'.format(cell.age), QTextOption(Qt.AlignCenter))
        
        QGraphicsPixmapItem.__init__(self, pixmap, *args, **kwargs)
        self.cell = cell
