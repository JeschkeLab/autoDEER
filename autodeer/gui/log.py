import typing
from PyQt6.QtWidgets import QApplication, QWidget,QLabel,QDoubleSpinBox,QGridLayout,QAbstractSpinBox,QDialog
from PyQt6 import QtCore, uic
import PyQt6.QtCore as QtCore
import PyQt6.QtGui as QtGui
import os
QtCore.QDir.addSearchPath('icons', 'gui/resources')
package_directory = os.path.dirname(os.path.abspath(__file__))


class LogDialog(QDialog):
    def __init__(self, parent: QWidget | None = ...) -> None:
        super().__init__(parent)

        uic.loadUi(os.path.join(package_directory, 'log.ui'), self)

        self.log_levels = {
            'DEBUG': 10,
            'INFO': 20,
            'WARNING': 30,
            'ERROR': 40,
            'CRITICAL': 50
        }

        self.autoDEER_log_level_combo.addItems(self.log_levels.keys())
        self.autoDEER_log_level_combo.setCurrentText('INFO')
        self.interface_log_level_combo.addItems(self.log_levels.keys())
        self.interface_log_level_combo.setCurrentText('INFO')

        self.autoDEER_log_level_combo.currentTextChanged.connect(self.update_log_level)
        self.interface_log_level_combo.currentTextChanged.connect(self.update_log_level)

        self.clear_log_button.clicked.connect(self.clear_log_text)


    def update_log_text():

        pass

    def clear_log_text():

        pass