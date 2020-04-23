from PySide2.QtCore import *
from PySide2.QtGui import *
from PySide2.QtWidgets import *

class Settings():
    show_protein_deletion_threshold = False
    show_cell_division_threshold = False
    show_sensor_reinforcement_threshold = False
    show_sym_prob_threshold = False

class SettingsArea(QWidget):
    settingChanged = Signal(str)
    
    def __init__(self, *args, **kwargs):
        QWidget.__init__(self, *args, **kwargs)

        layout = QFormLayout()

        protein_deletion_threshold_checkbox = QCheckBox()
        layout.addRow('Show Protein Deletion Threshold:', protein_deletion_threshold_checkbox)
        protein_deletion_threshold_checkbox.setCheckState(SettingsArea.bool_to_state(Settings.show_protein_deletion_threshold))
        protein_deletion_threshold_checkbox.stateChanged.connect(lambda state: self.update_bool_setting('show_protein_deletion_threshold', state))

        cell_division_threshold_checkbox = QCheckBox()
        layout.addRow('Show Cell Division Threshold:', cell_division_threshold_checkbox)
        cell_division_threshold_checkbox.setCheckState(SettingsArea.bool_to_state(Settings.show_cell_division_threshold))
        cell_division_threshold_checkbox.stateChanged.connect(lambda state: self.update_bool_setting('show_cell_division_threshold', state))

        sensor_reinforcement_threshold_checkbox = QCheckBox()
        layout.addRow('Show Sensor Reinforcement Threshold:', sensor_reinforcement_threshold_checkbox)
        sensor_reinforcement_threshold_checkbox.setCheckState(SettingsArea.bool_to_state(Settings.show_sensor_reinforcement_threshold))
        sensor_reinforcement_threshold_checkbox.stateChanged.connect(lambda state: self.update_bool_setting('show_sensor_reinforcement_threshold', state))
        

        sym_prob_threshold_checkbox = QCheckBox()
        layout.addRow('Show Sym Prob Threshold:', sym_prob_threshold_checkbox)
        sym_prob_threshold_checkbox.setCheckState(SettingsArea.bool_to_state(Settings.show_sym_prob_threshold))
        sym_prob_threshold_checkbox.stateChanged.connect(lambda state: self.update_bool_setting('show_sym_prob_threshold', state))

        self.setLayout(layout)

    @Slot()
    def update_bool_setting(self, name, check_state):
        setattr(Settings, name, SettingsArea.state_to_bool(check_state))
        self.settingChanged.emit(name)

    @staticmethod
    def state_to_bool(state):
        return state == Qt.Checked

    @staticmethod
    def bool_to_state(val):
        return Qt.Checked if val else Qt.Unchecked
