from PySide2.QtCore import *
from PySide2.QtGui import *
from PySide2.QtWidgets import *

class ToolbarArea(QToolBar):
    indexChanged = Signal(tuple)
    
    def __init__(self, run, *args, **kwargs):
        QToolBar.__init__(self, *args, **kwargs)
        
        self.eaStepSpin = QSpinBox()
        self.eaStepSpin.setRange(0, run.ea_steps)
        self.eaStepSpin.valueChanged.connect(self.handle_index_changed)
        #self.eaStepSpin.setSingleStep(run.step_interval.step)

        self.indivSpin = QSpinBox()
        self.indivSpin.setRange(1, run.pop_size)
        self.indivSpin.valueChanged.connect(self.handle_index_changed)
        
        self.regStepSpin = QSpinBox()
        self.regStepSpin.setRange(1, run.reg_steps + 1)
        self.regStepSpin.valueChanged.connect(self.handle_index_changed)

        self.addWidget(QLabel("EA Step:"))
        self.addWidget(self.eaStepSpin)
        self.addWidget(QLabel("Indiv:"))
        self.addWidget(self.indivSpin)
        self.addWidget(QLabel("Reg Step:"))
        self.addWidget(self.regStepSpin)

    def getIndex(self):
        return (self.eaStepSpin.value(), self.indivSpin.value(), self.regStepSpin.value())
        
    @Slot()
    def handle_index_changed(self, new_val):
        self.indexChanged.emit(self.getIndex())
    
