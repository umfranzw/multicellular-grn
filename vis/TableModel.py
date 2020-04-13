from PySide2.QtWidgets import *
from PySide2.QtGui import *
from PySide2.QtCore import *

class TableModel(QAbstractTableModel):
    rowChanged = Signal(str, bool)
    #note: a list of valid strings can be obtained by calling QColor.colorNames()
    colours = list(map(QColor, [
        'skyblue',
        'yellowgreen',
        'red',
        'darkviolet',
        'forestgreen',
        'deeppink',
        'orange',
        'navy',
        'slategray'
    ]))
    
    def __init__(self, data, headers):
        super(TableModel, self).__init__()
        self._data = data
        self._checks = {}
        self._headers = headers

    def checkState(self, index):
        if index in self._checks:
            return self._checks[index]
        else:
            return Qt.Unchecked

    def getCheckedInfo(self):
        #return a list of the form [(julia_props_obj, QColor), ...] that contains one entry for each checked row
        info = []
        for pindex in self._checks.keys():
            if self.checkState(pindex) == Qt.Checked:
                props = self._data[pindex.row()][-1]
                colour = self.get_colour(pindex.row())
                info.append((props, colour))

        return info

    def data(self, index, role=Qt.DisplayRole):
        if role == Qt.DisplayRole:
            return self._data[index.row()][index.column()]
        elif role == Qt.CheckStateRole:
            if index.column() == 0:
                #we use QPersistentModelIndex so that if the table changes, the ref still points to the same row
                return self.checkState(QPersistentModelIndex(index))
        elif role == Qt.DecorationRole:
            if index.column() == self.columnCount(index) - 1:
                colour = self.get_colour(index.row())
                pixmap = QPixmap(30, 20)
                pixmap.fill(colour)
                return pixmap
        # elif role == Qt.BackgroundColorRole:
        #     if index.column() == self.columnCount(index) - 1: 
        #         return TableModel.colours[index.row() % len(TableModel.colours)]
    
        return None

    def get_colour(self, row):
        return TableModel.colours[row % len(TableModel.colours)]

    def setData(self, index, value, role):
        if not index.isValid():
            return False
        if role == Qt.CheckStateRole:
            self._checks[QPersistentModelIndex(index)] = value
            self.rowChanged.emit(self._data[index.row()][index.column()], value)
            return True
        return False
            
    def rowCount(self, index):
        return len(self._data)

    def columnCount(self, index):
        return len(self._data[0]) #(all rows are of equal length, +1 for colour column, -1 for pointer at the end)

    def headerData(self, section, orientation, role):
        if role == Qt.DisplayRole:
            if orientation == Qt.Horizontal:
                return self._headers[section]
            else:
                return None

    def flags(self, index):
        flag = QAbstractTableModel.flags(self, index)
        if index.column() == 0:
            flag |= Qt.ItemIsEditable | Qt.ItemIsUserCheckable
            
        return flag
