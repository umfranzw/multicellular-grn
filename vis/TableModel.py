from PySide2.QtWidgets import *
from PySide2.QtGui import *
from PySide2.QtCore import *
import colorsys

class TableModel(QAbstractTableModel):
    rowChanged = Signal(list)

    kelly_colour_vals = (0x222222, 0xf3c300, 0x875692, 0xf38400, 0xa1caf1, 0xbe0032, 0xc2b280, 0x848482, 0x008856, 0xe68fac, 0x0067a5, 0xf99379, 0x604e97, 0xf6a600, 0xb3446c, 0xdcd300, 0x882d17, 0x8db600, 0x654522, 0xe25822, 0x2b3d26) #note: white has been removed

    colours = [QColor.fromRgb(c) for c in kelly_colour_vals]
    
    def __init__(self, data, headers):
        super(TableModel, self).__init__()
        self._data = data
        self._checks = {} # julia hash value (from _data[-1] -> (Qt.Checked/UnChecked, julia props obj)
        self._headers = headers

    def refresh(self, data):
        new_checks = {}
        for key in self._checks.keys():
            if self.checkState(key) == Qt.Checked:
                new_checks[key] = self._checks[key]

        self.beginResetModel()
        self._checks = new_checks
        self._data = data
        self.endResetModel()

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
                colour = self.getColour(key)
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
                colour = self.getColour(self._data[index.row()][-1])
                pixmap = QPixmap(30, 20)
                pixmap.fill(colour)
                return pixmap
    
        return None

    def getColour(self, val):
        colour_index = val % len(TableModel.colours)
        return TableModel.colours[colour_index]

    def setData(self, index, value, role):
        if not index.isValid():
            return False
        if role == Qt.CheckStateRole:
            self._checks[self._data[index.row()][-1]] = (value, self._data[index.row()][-2])
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
