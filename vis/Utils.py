from PySide2.QtGui import *
from PySide2.QtWidgets import *
from PySide2.QtCore import *
from PySide2.QtCharts import *

class Utils():
    @staticmethod
    @Slot()
    def save_chart_view(chart_view):
        filename = Utils.run_save_image_dialog()
        pixmap = QPixmap(chart_view.size())
        chart_view.render(pixmap)
        pixmap.save(filename, format='png')

    @Slot()
    @staticmethod
    def save_graphics_view(graphics_view):
        filename = Utils.run_save_image_dialog()
        pixmap = QPixmap(graphics_view.size())
        graphics_view.render(pixmap)
        pixmap.save(filename, format='png')

    @Slot()
    @staticmethod
    def save_pixmap(pixmap):
        filename = Utils.run_save_image_dialog()
        pixmap.save(filename, format='png')


    @staticmethod
    def run_save_image_dialog():
        filename = QFileDialog.getSaveFileName(
            caption='Save File',
            filter='PNG Image (*.png)'
        )
        
        filename = filename[0]
        if not filename.endswith('.png'):
            filename += '.png'
            
        return filename
