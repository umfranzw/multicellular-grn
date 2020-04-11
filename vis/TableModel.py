from PySide2.QtWidgets import *
from PySide2.QtGui import *
from PySide2.QtCore import *

class TableModel(QAbstractTableModel):
    rowChanged = Signal(str, bool)
    
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

    def getCheckedProps(self):
        checked = []
        for pindex in self._checks.keys():
            if self.checkState(pindex) == Qt.Checked:
                checked.append(self._data[pindex.row()][-1])

        return checked

    def data(self, index, role=Qt.DisplayRole):
        if role == Qt.DisplayRole:
            return self._data[index.row()][index.column()]
        elif role == Qt.CheckStateRole:
            if index.column() == 0:
                #we use QPersistentModelIndex so that if the table changes, the ref still points to the same row
                return self.checkState(QPersistentModelIndex(index))
        return None

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
        return len(self._data[0]) - 1 #(all rows are of equal length, -1 for pointer at the end)

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
