from PySide2.QtCore import *
from PySide2.QtGui import *
from PySide2.QtWidgets import *

#https://stackoverflow.com/questions/40746350/why-qspinbox-jumps-twice-the-step-value
class CustomStyle(QProxyStyle):
    def styleHint(self, hint, option=None, widget=None, returnData=None):
        if hint == QStyle.SH_SpinBox_ClickAutoRepeatThreshold:
            return 2**31 - 1 #max 32-bit int value (since this'll be passed to a C++ library)
        else:
            return super().styleHint(hint, option, widget, returnData)
