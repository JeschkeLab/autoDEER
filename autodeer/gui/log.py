import typing
from PyQt6.QtWidgets import QApplication, QWidget,QLabel,QFileDialog,QGridLayout,QAbstractSpinBox,QDialog
from PyQt6 import QtCore, uic
import PyQt6.QtCore as QtCore
import PyQt6.QtGui as QtGui
import os
import logging
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

        self.autoDEER_log_level_combo.currentTextChanged.connect(self.update_log_text)
        self.interface_log_level_combo.currentTextChanged.connect(self.update_log_text)

        self.clear_log_button.clicked.connect(self.clear_log)
        self.close_log_button.clicked.connect(self.hide)
        self.save_log_button.clicked.connect(self.save_log_to_file)


        core_logger = logging.getLogger('autoDEER')
        hardware_logger = logging.getLogger('interface')

        core_logger.handlers[1].signal.connect(self.add_new_entry)

        self.log = []

    def add_new_entry(self, entry: dict):
        self.log.append(entry)

        if len(self.log) > 1000:
            self.log.pop(0)

        self.update_log_text()

    def update_log_text(self):
        
        formatter = "%(time)s [%(name)s] - %(level)s: %(message)s \n"
        self.log_text.clear()
        text = ''

        autoDEER_filter_level = self.log_levels[self.autoDEER_log_level_combo.currentText()]
        interface_filter_level = self.log_levels[self.interface_log_level_combo.currentText()]

        for entry in self.log:

            if ('autoDEER' in entry['name']) and (self.log_levels[entry['level']] >= autoDEER_filter_level):
                text += (formatter % entry)
            elif ('interface' in entry['name']) and (self.log_levels[entry['level']] >= interface_filter_level):
                text += (formatter % entry)

        self.log_text.insertPlainText(text)
        pass
    def clear_log(self):
        self.log = []
        self.update_log_text()
        pass

    def save_log_to_file(self):
        save_path = QFileDialog.getSaveFileName(self, 'Save File', self.current_folder, ".txt")

        if save_path[0] != '':
            with open(save_path[0], 'w') as f:
                f.write(self.log_text.toPlainText())
        pass