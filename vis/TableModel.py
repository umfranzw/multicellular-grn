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

    def refresh(self, data):
        # checked = set()
        # for pindex in self._checks.keys():
        #     if self.checkState(pindex) == Qt.Checked:
        #         props = tuple(self._data[pindex.row()][:-1]) #get all properties except the julia pointer on the end
        #         checked.add(props)

        self.beginResetModel()
        self._checks = {}

        # for row in range(len(data)):
        #     key = tuple(data[row][:-1])
        #     if key in checked:
        #         index = self.index(row, 0)
        #         self._checks[QPersistentModelIndex(index)] = Qt.Checked

        #self.removeRows(0, self.rowCount())

        self._data = data
        self.endResetModel()

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
                colour = self.getColour(pindex.row())
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
                colour = self.getColour(index.row())
                pixmap = QPixmap(30, 20)
                pixmap.fill(colour)
                return pixmap
    
        return None

    def getColour(self, row):
        key = tuple(self._data[row][:-1]) #get everything except the julia pointer on the end
        colour_index = hash(key) % len(TableModel.colours)
        return TableModel.colours[colour_index]

    def setData(self, index, value, role):
        if not index.isValid():
            return False
        if role == Qt.CheckStateRole:
            self._checks[QPersistentModelIndex(index)] = value
            self.rowChanged.emit(self._data[index.row()][index.column()], value)
            return True
        return False
            
    def rowCount(self, index=QModelIndex()):
        return len(self._data)

    def columnCount(self, index=QModelIndex()):
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
