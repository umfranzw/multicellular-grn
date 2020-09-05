from PySide2.QtCore import *
from PySide2.QtGui import *
from PySide2.QtWidgets import *

from CustomStyle import CustomStyle

class ToolbarArea(QToolBar):
    indexChanged = Signal(tuple) #fires when anything changes
    eaStepChanged = Signal(int)
    showBestChanged = Signal(bool)
    
    def __init__(self, run, data_tools, *args, **kwargs):
        QToolBar.__init__(self, *args, **kwargs)
        
        self.eaStepSpin = QSpinBox()
        self.eaStepSpin.setStyle(CustomStyle())
        self.eaStepSpin.setRange(run.step_range.start, run.step_range.stop - 1) #note: in Julia the top bound is inclusive, in Python it's exclusive
        self.eaStepSpin.setSingleStep(run.step_range.step)
        self.eaStepSpin.valueChanged.connect(self.handle_index_changed)
        self.eaStepSpin.valueChanged.connect(self.handle_ea_step_changed)

        self.indivSpin = QSpinBox()
        self.indivSpin.setStyle(CustomStyle())
        self.indivSpin.setRange(1, run.pop_size)
        self.indivSpin.valueChanged.connect(self.handle_index_changed)
        
        self.regStepSpin = QSpinBox()
        self.regStepSpin.setStyle(CustomStyle())
        self.regStepSpin.setRange(1, run.reg_steps + 1)
        self.regStepSpin.valueChanged.connect(self.handle_index_changed)

        self.show_best_checkbox = QCheckBox()
        self.show_best_checkbox.stateChanged.connect(self.handle_show_best_changed)

        self.seed_label = QLabel('Seed: {}'.format(data_tools.get_base_seed()))
        self.seed_label.setTextInteractionFlags(Qt.TextSelectableByMouse)

        self.addWidget(QLabel("EA Step:"))
        self.addWidget(self.eaStepSpin)
        self.addWidget(QLabel("Indiv:"))
        self.addWidget(self.indivSpin)
        self.addWidget(QLabel("Reg Step:"))
        self.addWidget(self.regStepSpin)
        self.addWidget(QLabel("Show Best Captured"))
        self.addWidget(self.show_best_checkbox)
        self.addWidget(self.seed_label)

    def getIndex(self):
        return (self.eaStepSpin.value(), self.indivSpin.value(), self.regStepSpin.value())

    def setIndex(self, index):
        ea_step, pop_index, reg_step = index

        #self.blockSignals(True)
        #ea_step_changed = self.eaStepSpin.value() != ea_step
        self.eaStepSpin.setValue(ea_step)
        #index_changed = ea_step_changed or self.indivSpin.value() != pop_index
        self.indivSpin.setValue(pop_index)
        #index_changed = index_changed or self.regStepSpin.value() != reg_step
        self.regStepSpin.setValue(reg_step)
        #self.blockSignals(False)

        #fire each signal *once* (if necessary) to update the view
        # if index_changed:
        #     self.handle_index_changed(reg_step) #assume they've just changed the reg_step (not great but works)
        # if ea_step_changed:
        #     self.handle_ea_step_changed(ea_step)

    def disable_spin_buttons(self):
        self.set_enabled_state(False)

    def enable_spin_buttons(self):
        self.set_enabled_state(True)

    def set_enabled_state(self, state):
        self.eaStepSpin.setEnabled(state)
        self.indivSpin.setEnabled(state)
        self.regStepSpin.setEnabled(state)
    
    @Slot()
    def handle_index_changed(self, new_val):
        self.indexChanged.emit(self.getIndex())

    @Slot()
    def handle_ea_step_changed(self, new_val):
        self.eaStepChanged.emit(new_val)

    @Slot()
    def handle_show_best_changed(self, new_val):
        self.showBestChanged.emit(new_val == Qt.Checked)
    
