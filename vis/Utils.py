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
        pixmap = QPixmap(graphics_view.scene.width(), graphics_view.scene.height())
        pixmap.fill(Qt.white)
        painter = QPainter(pixmap)
        graphics_view.scene.render(painter)
        painter.end()
        pixmap.save(filename, format='png')

    @Slot()
    @staticmethod
    def save_pixmap(pixmap):
        filename = Utils.run_save_image_dialog()
        pixmap.save(filename, format='png')


    @staticmethod
    def run_save_image_dialog():
        return Utils.run_save_file_dialog('PNG Image (*.png)', '.png')

    @staticmethod
    def run_save_csv_dialog():
        return Utils.run_save_file_dialog('Comma Separated Value (*.csv)', '.csv')

    @staticmethod
    def run_save_file_dialog(filter_text, file_suffix):
        filename = QFileDialog.getSaveFileName(
            caption='Save File',
            filter=filter_text
        )
        
        filename = filename[0]
        if not filename.endswith(file_suffix):
            filename += file_suffix
            
        return filename

    @staticmethod
    def run_save_folder_dialog():
        path = QFileDialog.getExistingDirectory(
            caption='Save All Steps'
        )
        
        return path

    
