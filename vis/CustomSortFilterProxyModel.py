from PySide2.QtCore import QSortFilterProxyModel, QRegExp, Qt

class CustomSortFilterProxyModel(QSortFilterProxyModel):
    def __init__(self, parent):
        QSortFilterProxyModel.__init__(self, parent)

    def filterAcceptsRow(self, sourceRow, sourceParent):
        sourceModel = self.sourceModel()
        filter_str = self.filterRegExp().pattern()
        components = filter_str.split(',')
        
        #any unspecified columns at the end are assumed to match
        accept = True
        col = 0
        while accept and col < len(components):
            index = sourceModel.index(sourceRow, col, sourceParent)
            regex = QRegExp('^' + components[col], cs=Qt.CaseSensitivity.CaseInsensitive)
            
            #convert to a string to simplify regexp matching
            data_item = str(sourceModel.data(index))
            accept = accept and regex.indexIn(data_item) != -1
            col += 1

        return accept
        
