from PySide2.QtWidgets import *
from PySide2.QtGui import *
from PySide2.QtCore import *

from ColourPicker import ColourPicker

class TableModel(QAbstractTableModel):
    rowChanged = Signal(list)

    def __init__(self, data, headers):
        super(TableModel, self).__init__()
        self._data = data
        self._checks = {} # julia hash value (from _data[row][-1] -> (Qt.Checked/UnChecked, julia props obj)
        self._headers = headers
        self.colour_picker = ColourPicker()

    def reset_colour_picker(self):
        self.colour_picker.reset()
        
    def refresh(self, data):
        self.beginResetModel()
        self._data = data
        self.endResetModel()

    def setCheckedRows(self, row_indices):
        self.beginResetModel()
        self._checks.clear()
        for row in row_indices:
            key = self._data[row][-1]
            value = (Qt.Checked, self._data[row][-2])
            self._checks[key] = value
        self.endResetModel()
        
        self.rowChanged.emit(self.getCheckedInfo())

    def checkState(self, index):
        if index in self._checks:
            return self._checks[index][0]
        else:
            return Qt.Unchecked

    def getCheckedInfo(self):
        #return a list of the form [(julia_props_obj, QColor), ...] that contains one entry for each checked row
        info = []
        for key in self._checks.keys():
            if self.checkState(key) == Qt.Checked:
                props = self._checks[key][1]
                colour = self.colour_picker.get(key)
                info.append((props, colour))

        return info

    def data(self, index, role=Qt.DisplayRole):
        if role == Qt.DisplayRole:
            return self._data[index.row()][index.column()]
        elif role == Qt.CheckStateRole:
            if index.column() == 0:
                return self.checkState(self._data[index.row()][-1])
        elif role == Qt.DecorationRole:
            if index.column() == self.columnCount(index) - 1:
                colour = self.colour_picker.get(self._data[index.row()][-1])
                pixmap = QPixmap(30, 20)
                pixmap.fill(colour)
                return pixmap
    
        return None

    def setData(self, index, value, role, emit=True):
        if not index.isValid():
            return False
        if role == Qt.CheckStateRole:
            self._checks[self._data[index.row()][-1]] = (value, self._data[index.row()][-2])
            if emit:
                self.rowChanged.emit(self.getCheckedInfo())
            return True
        
        return False
            
    def rowCount(self, index=QModelIndex()):
        return len(self._data)

    def columnCount(self, index=QModelIndex()):
        return len(self._headers)

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
