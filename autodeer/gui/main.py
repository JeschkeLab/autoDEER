from PyQt6.QtWidgets import QApplication, QMainWindow, QFileDialog,QMessageBox, QDialog, QPushButton,QVBoxLayout, QLineEdit,QHBoxLayout, QLabel, QProgressDialog
import PyQt6.QtWidgets as QtWidgets
from PyQt6 import uic
import PyQt6.QtCore as QtCore
import PyQt6.QtGui as QtGui
from PyQt6.QtSvgWidgets import QSvgWidget
from qt_material import apply_stylesheet
from pathlib import Path
import sys, traceback, os
from threadpoolctl import threadpool_limits


from matplotlib.backends.backend_qtagg import FigureCanvas, NavigationToolbar2QT
from matplotlib.gridspec import GridSpec

import matplotlib.pyplot as plt
import autodeer as ad
import pyepr as epr
import numpy as np
from autodeer.gui.tools import *
from autodeer.gui.autoDEER_worker import autoDEERWorker
from autodeer.gui.deer_panel import DEERplot
from autodeer.gui.log import LogDialog
from autodeer.gui.modetune import ModeTune
from autodeer.gui.constants import *
import yaml
import time
import datetime
import logging
from autodeer.Logging import setup_logs, change_log_level
from deerlab import store_pickle
import deerlab as dl
import copy
main_log = logging.getLogger('autoDEER')
from queue import Queue

package_directory = os.path.dirname(os.path.abspath(__file__))

QtCore.QDir.addSearchPath('icons', f"{package_directory}/resources")

BackgroundModels = {'hom3D':dl.bg_hom3d, 'hom3d ex':dl.bg_hom3dex, 'hom3d frac':dl.bg_homfractal,'str exp':dl.bg_strexp}

class WorkerSignals(QtCore.QObject):
    '''
    Defines the signals available from a running worker thread.

    Supported signals are:

    finished
        No data

    error
        tuple (exctype, value, traceback.format_exc() )

    result
        object data returned from processing, anything

    progress
        int indicating % progress

    '''
    finished = QtCore.pyqtSignal()
    error = QtCore.pyqtSignal(tuple)
    result = QtCore.pyqtSignal(object)
    progress = QtCore.pyqtSignal(int)

class Worker(QtCore.QRunnable):
    '''
    Worker thread

    Inherits from QRunnable to handler worker thread setup, signals and wrap-up.

    :param callback: The function callback to run on this worker thread. Supplied args and
                     kwargs will be passed through to the runner.
    :type callback: function
    :param args: Arguments to pass to the callback function
    :param kwargs: Keywords to pass to the callback function

    '''

    def __init__(self, fn, *args, **kwargs):
        super(Worker, self).__init__()

        # Store constructor arguments (re-used for processing)
        self.fn = fn
        self.args = args
        self.kwargs = kwargs
        self.signals = WorkerSignals()

        # # Add the callback to our kwargs
        # self.kwargs['progress_callback'] = self.signals.progress

    @QtCore.pyqtSlot()
    def run(self):
        '''
        Initialise the runner function with passed args, kwargs.
        '''

        # Retrieve args/kwargs here; and fire processing using them
        try:
            result = self.fn(*self.args, **self.kwargs)
        except:
            traceback.print_exc()
            exctype, value = sys.exc_info()[:2]
            self.signals.error.emit((exctype, value, traceback.format_exc()))
        else:
            self.signals.result.emit(result)  # Return the result of the processing
        finally:
            self.signals.finished.emit()  # Done

def fieldsweep_fit(fsweep_analysis,fit):
    if fit:
        fsweep_analysis.fit(xtol=1e-5, lin_maxiter=100)
    fsweep_analysis.smooth() # Always smooth anyway
    return fsweep_analysis

def respro_process(dataset,f_lims, fieldsweep=None,cores=1):
    respro = epr.ResonatorProfileAnalysis(
        dataset,f_lims=f_lims
    )
    fc_guess = respro.freqs[respro.freqs.shape[0]//2]
    with threadpool_limits(limits=cores, user_api='blas'):
        respro.fit(cores=cores,fc_guess=fc_guess)

    if fieldsweep is not None:
        LO_new = fieldsweep.freq + epr.optimise_spectra_position(respro, fieldsweep)
        return respro, LO_new


    return respro

def relax_process(dataset):

    seq_name = dataset.seq_name

    if (seq_name == 'CarrPurcellSequence') :
        if dataset['tau1'].max() > 1e4:
            dataset['tau1'] /= 1e3
        CP_data = epr.CarrPurcellAnalysis(dataset)
        CP_data.fit('auto')
        return CP_data
    
    elif seq_name == 'T2RelaxationSequence':
        Tm_data = epr.HahnEchoRelaxationAnalysis(dataset)
        Tm_data.fit('auto')

        return Tm_data
    
    elif seq_name == 'RefocusedEcho1DSequence':
        Tm_data = ad.RefocusedEcho1DAnalysis(dataset)
        Tm_data.fit('auto')

        return Tm_data


    

def T2_process(dataset):
    Tm_data = epr.HahnEchoRelaxationAnalysis(dataset)
    Tm_data.fit('auto')

    return Tm_data

class autoDEERUI(QMainWindow):


    def __init__(self):
        super().__init__()
 
        # loading the ui file with uic module
        uic.loadUi(f"{package_directory}/gui2.ui", self)
        logo_pixmap = QtGui.QPixmap('icons:logo.png')
        logo_pixmap.setDevicePixelRatio(1.0)
        logo_pixmap = logo_pixmap.scaledToHeight(60)
        self.logo.setPixmap(logo_pixmap)
        self.set_spectrometer_connected_light(0)

        self.setWindowTitle("autoDEER")
        self.Version_label.setText(f"Version: {ad.__version__}")


        self.qDEER_tab.layout().addWidget(DEERplot())
        self.q_DEER = self.qDEER_tab.layout().itemAt(0).widget()
        self.lDEER_tab.layout().addWidget(DEERplot())
        self.longDEER = self.lDEER_tab.layout().itemAt(0).widget()


        self.create_setup_figure()
        self.create_relax_figure()
        # self.advanced_mode_inputs()
        self.build_advanced_mode()

        self.threadpool = QtCore.QThreadPool()
        self.current_results = {}
        self.current_data = {}
        self.pulses = {}
        self.spectromterInterface = None
        self.waitCondition = None
        self.queue = Queue()

        self.FullyAutoButton.clicked.connect(self.RunFullyAutoDEER)
        self.AdvancedAutoButton.clicked.connect(self.RunAdvancedAutoDEER)

        docs_url = QtCore.QUrl('https://jeschkelab.github.io/autoDEER/')
        github_url = QtCore.QUrl('https://github.com/JeschkeLab/autoDEER/')
        issues_url = QtCore.QUrl('https://github.com/JeschkeLab/autoDEER/issues')
        discussion_url = QtCore.QUrl('https://github.com/JeschkeLab/autoDEER/discussions')

        self.actionDocumentation.triggered.connect(lambda: QtGui.QDesktopServices.openUrl(docs_url))
        self.actionGitHub.triggered.connect(lambda: QtGui.QDesktopServices.openUrl(github_url))
        self.actionIssues.triggered.connect(lambda: QtGui.QDesktopServices.openUrl(issues_url))
        self.actionDiscussions.triggered.connect(lambda: QtGui.QDesktopServices.openUrl(discussion_url))
        self.actionOpen_Folder.triggered.connect(self.load_folder)

        self.actionLoadConfig.triggered.connect(lambda : self.load_spectrometer_config(None))
        self.actionConnect.triggered.connect(self.connect_spectrometer)
        self.actionSaveReport.triggered.connect(self.create_report)
        
        self.show_respro.clicked.connect(lambda: self.resonatorProfileFigure())
        self.show_EDFS.clicked.connect(lambda: self.refresh_fieldsweep_after_fit())
        self.OptimisePulsesButton.clicked.connect(lambda: self.optimise_pulses_button())
        

        self.Tab_widget.setCurrentIndex(0)
        # set current folder to home directory
        self.current_folder = None
        self.config = None
        self.connected = False
        self.DL_params = {}
        self.worker = None
        self.starttime = time.time()
        self.Bruker=False
        self.userinput = {'label_eff':100,'MaxTime':24}

        self.LO = 0
        self.Q = 0
        self.gyro = 0.002803632236095
        self.gyro = 0.002808859721083
        self.cores = 1
        self.Min_tp=12

        self.deer_settings = {'ESEEM':None, 'ExpType':'5pDEER'}
        self.priorties = {'Auto': 150, 'MNR':300, 'Distance': 80, 'Single':200}

        self.priotityComboBox.addItems(list(self.priorties.keys()))
        self.correction_factor=1
        self.est_lambda = None
        self.pump_pulses = []

    def build_advanced_mode(self):
        
        # Link Checkbox to enbaling group boxes
        self.AdvSequenceCheck.stateChanged.connect(lambda state: self.SequenceGroupBox.setEnabled(state == QtCore.Qt.CheckState.Checked.value))
        self.AdvPulsesCheck.stateChanged.connect(lambda state: self.PulsesGroupBox.setEnabled(state == QtCore.Qt.CheckState.Checked.value))
        self.AdvSequenceCheck.setChecked(False)
        self.AdvPulsesCheck.setChecked(False)
        self.PulsesGroupBox.setEnabled(False)
        self.SequenceGroupBox.setEnabled(False)

        self.AdvModeSave.clicked.connect(self.save_advaced_options)
        self.AdvModeLoad.clicked.connect(self.load_advanced_options)

        # Populate Sequence Options
        self.Exp_types.addItems(['auto', '5pDEER', '4pDEER', ]) # 'Ref2D'
        self.bg_model_combo.addItems(list(BackgroundModels.keys()))
        self.AdvSeqOptions = {}
        self.Exp_types.currentIndexChanged.connect(self.update_advanced_sequence_options)
        self.update_advanced_sequence_options()


        # Populate Pulse Options Table
        pulse_table:QtWidgets.QTableWidget  = self.PulseParamTable
        # Make the columns expand to fill the space
        pulse_table.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.ResizeMode.Stretch)

        self.AdvPulseSelect = {}
        # Row 0: Headers
        # Row 1: Type

        self.AdvPulseSelect['ExcType'] = QtWidgets.QComboBox()
        self.AdvPulseSelect['ExcType'].addItems(['Auto', 'Rect.', 'Gauss'])
        pulse_table.setCellWidget(1,1,self.AdvPulseSelect['ExcType'])

        self.AdvPulseSelect['RefType'] = QtWidgets.QComboBox()
        self.AdvPulseSelect['RefType'].addItems(['Auto', 'Rect.', 'Gauss'])
        pulse_table.setCellWidget(1,2,self.AdvPulseSelect['RefType'])
        
        self.AdvPulseSelect['PumpType'] = QtWidgets.QComboBox()
        self.AdvPulseSelect['PumpType'].addItems(['Auto', 'Rect.', 'Chirp', 'HS', 'Gauss'])
        pulse_table.setCellWidget(1,3,self.AdvPulseSelect['PumpType'])

        # Row 2: Center Freq
        self.AdvPulseSelect['ExcFreq'] = QtWidgets.QDoubleSpinBox()
        self.AdvPulseSelect['ExcFreq'].setRange(-1000,1000)
        self.AdvPulseSelect['ExcFreq'].setSuffix(' MHz')
        pulse_table.setCellWidget(2,1,self.AdvPulseSelect['ExcFreq'])
        
        self.AdvPulseSelect['RefFreq'] = QtWidgets.QDoubleSpinBox()
        self.AdvPulseSelect['RefFreq'].setRange(-1000,1000)
        self.AdvPulseSelect['RefFreq'].setSuffix(' MHz')
        pulse_table.setCellWidget(2,2,self.AdvPulseSelect['RefFreq'])
        
        self.AdvPulseSelect['PumpFreq'] = QtWidgets.QDoubleSpinBox()
        self.AdvPulseSelect['PumpFreq'].setRange(-1000,1000)
        self.AdvPulseSelect['PumpFreq'].setSuffix(' MHz')
        pulse_table.setCellWidget(2,3,self.AdvPulseSelect['PumpFreq'])
        # Row 3. Pulse Length
        self.AdvPulseSelect['ExcLength'] = QtWidgets.QDoubleSpinBox()
        self.AdvPulseSelect['ExcLength'].setRange(0,500)
        self.AdvPulseSelect['ExcLength'].setSuffix(' ns')
        pulse_table.setCellWidget(3,1,self.AdvPulseSelect['ExcLength'])
        
        self.AdvPulseSelect['RefLength'] = QtWidgets.QDoubleSpinBox()
        self.AdvPulseSelect['RefLength'].setRange(0,500)
        self.AdvPulseSelect['RefLength'].setSuffix(' ns')
        pulse_table.setCellWidget(3,2,self.AdvPulseSelect['RefLength'])
        self.AdvPulseSelect['PumpLength'] = QtWidgets.QDoubleSpinBox()
        self.AdvPulseSelect['PumpLength'].setRange(0,500)
        self.AdvPulseSelect['PumpLength'].setSuffix(' ns')
        pulse_table.setCellWidget(3,3,self.AdvPulseSelect['PumpLength'])
        
        # Row 4. Bandwidth
        self.AdvPulseSelect['ExcBandwidth'] = QtWidgets.QDoubleSpinBox()
        self.AdvPulseSelect['ExcBandwidth'].setRange(0,1000)
        self.AdvPulseSelect['ExcBandwidth'].setSuffix(' MHz')
        pulse_table.setCellWidget(4,1,self.AdvPulseSelect['ExcBandwidth'])
        
        self.AdvPulseSelect['RefBandwidth'] = QtWidgets.QDoubleSpinBox()
        self.AdvPulseSelect['RefBandwidth'].setRange(0,1000)
        self.AdvPulseSelect['RefBandwidth'].setSuffix(' MHz')
        pulse_table.setCellWidget(4,2,self.AdvPulseSelect['RefBandwidth'])
        
        self.AdvPulseSelect['PumpBandwidth'] = QtWidgets.QDoubleSpinBox()
        self.AdvPulseSelect['PumpBandwidth'].setRange(0,1000)
        self.AdvPulseSelect['PumpBandwidth'].setSuffix(' MHz')
        pulse_table.setCellWidget(4,3,self.AdvPulseSelect['PumpBandwidth'])

        # if type if 'Auto' disable other options when changed
        self.AdvPulseSelect['ExcType'].currentIndexChanged.connect(lambda: self.update_advanced_pulse_options('Exc'))
        self.AdvPulseSelect['RefType'].currentIndexChanged.connect(lambda: self.update_advanced_pulse_options('Ref'))
        self.AdvPulseSelect['PumpType'].currentIndexChanged.connect(lambda: self.update_advanced_pulse_options('Pump'))
        self.AdvPulseSelect['ExcType'].setCurrentIndex(0)
        self.AdvPulseSelect['RefType'].setCurrentIndex(0)
        self.AdvPulseSelect['PumpType'].setCurrentIndex(0)
        self.update_advanced_pulse_options('Exc')
        self.update_advanced_pulse_options('Ref')
        self.update_advanced_pulse_options('Pump')

    def update_advanced_pulse_options(self, pulse_key):
        pulse_type = self.AdvPulseSelect[f'{pulse_key}Type'].currentText()
        if pulse_type == 'Auto':
            self.AdvPulseSelect[f'{pulse_key}Freq'].setEnabled(False)
            self.AdvPulseSelect[f'{pulse_key}Length'].setEnabled(False)
            self.AdvPulseSelect[f'{pulse_key}Bandwidth'].setEnabled(False)
        else:
            self.AdvPulseSelect[f'{pulse_key}Freq'].setEnabled(True)
            self.AdvPulseSelect[f'{pulse_key}Length'].setEnabled(True)
            if pulse_type in ['Chirp','HS']:
                self.AdvPulseSelect[f'{pulse_key}Bandwidth'].setEnabled(True)
            else:
                self.AdvPulseSelect[f'{pulse_key}Bandwidth'].setEnabled(False)

    def update_advanced_sequence_options(self):
        
        exp_type = self.Exp_types.currentText()
        self.SequenceSVGWidget:QSvgWidget
        # Clear previous widgets in SequenceParamRow
        while self.SequenceParamRow.count():
            item = self.SequenceParamRow.takeAt(0)
            widget = item.widget()
            if widget is not None:
                widget.deleteLater()

        if exp_type.lower() == 'auto':
            self.SequenceSVGWidget.setEnabled(False)
            # Clear SVG
            self.SequenceSVGWidget.hide()
           
        elif exp_type.lower() == '4pdeer':
            self.SequenceSVGWidget.setEnabled(True)
            self.SequenceSVGWidget.load('icons:4pDEERSeq.svg')
            self.SequenceSVGWidget.renderer().setAspectRatioMode(QtCore.Qt.AspectRatioMode.KeepAspectRatio)
            self.SequenceSVGWidget.show()
            for w in ['tau1','tau2','dt']:
                # Add QSpinBox to SequenceParamRow
                spinbox = QtWidgets.QDoubleSpinBox()
                spinbox.setRange(0,1000)
                if w == 'dt':
                    spinbox.setSuffix(' ns')
                    spinbox.setDecimals(1)
                    spinbox.setSingleStep(0.5)
                else:
                    spinbox.setSuffix(' us')
                    spinbox.setDecimals(2)
                    spinbox.setSingleStep(0.25)
                if w in self.AdvSeqOptions: # Preserve previous value if exists
                    spinbox.setValue(self.AdvSeqOptions[w])
                label = QtWidgets.QLabel(w)
                self.SequenceParamRow.addWidget(label)
                self.SequenceParamRow.addWidget(spinbox)
        elif exp_type.lower() == '5pdeer':
            self.SequenceSVGWidget.setEnabled(True)
            svg_path = QtCore.QDir.searchPaths('icons')[0] + '/5pDEERSeq.svg'
            self.SequenceSVGWidget.load('icons:5pDEERSeq.svg')
            self.SequenceSVGWidget.renderer().setAspectRatioMode(QtCore.Qt.AspectRatioMode.KeepAspectRatio)
            self.SequenceSVGWidget.show()
            for w in ['tau1','tau2','tau3','dt']:
                # Add QSpinBox to SequenceParamRow
                spinbox = QtWidgets.QDoubleSpinBox()
                spinbox.setRange(0,1000)
                if w == 'dt':
                    spinbox.setSuffix(' ns')
                    spinbox.setDecimals(1)
                    spinbox.setSingleStep(0.5)
                else:
                    spinbox.setSuffix(' us')
                    spinbox.setDecimals(2)
                    spinbox.setSingleStep(0.25)
                if w in self.AdvSeqOptions: # Preserve previous value if exists
                    spinbox.setValue(self.AdvSeqOptions[w])
                label = QtWidgets.QLabel(w)
                self.SequenceParamRow.addWidget(label)
                self.SequenceParamRow.addWidget(spinbox)
        else:
            raise NotImplementedError(f"Advanced sequence options for {exp_type} not implemented yet.")

        self.read_advanced_sequence_options()

    def read_advanced_sequence_options(self):
        exp_type = self.Exp_types.currentText()
        self.AdvSeqOptions = {}
        self.AdvSeqOptions['ExpType'] = exp_type
        if exp_type.lower() == 'auto':
            return
        elif exp_type == '4pDEER':
            self.AdvSeqOptions['tau1'] = self.SequenceParamRow.itemAt(1).widget().value()
            self.AdvSeqOptions['tau2'] = self.SequenceParamRow.itemAt(3).widget().value()
            self.AdvSeqOptions['dt'] = self.SequenceParamRow.itemAt(5).widget().value()
        elif exp_type == '5pDEER':
            self.AdvSeqOptions['tau1'] = self.SequenceParamRow.itemAt(1).widget().value()
            self.AdvSeqOptions['tau2'] = self.SequenceParamRow.itemAt(3).widget().value()
            self.AdvSeqOptions['tau3'] = self.SequenceParamRow.itemAt(5).widget().value()
            self.AdvSeqOptions['dt'] = self.SequenceParamRow.itemAt(7).widget().value()

    def read_advanced_pulse_options(self):
        self.AdvPulseOptions = {}
        for key, widget in self.AdvPulseSelect.items():
            if isinstance(widget, QtWidgets.QComboBox):
                self.AdvPulseOptions[key] = widget.currentText()
            elif isinstance(widget, QtWidgets.QDoubleSpinBox):
                if widget.suffix() == ' ns':
                    self.AdvPulseOptions[key] = widget.value()
                elif widget.suffix() == ' MHz': # Convert to GHz
                    self.AdvPulseOptions[key] = widget.value() / 1000.0
                else:
                    self.AdvPulseOptions[key] = widget.value()
    
    def validate_advanced_pulse_options(self):
        """
        Checks that if a pulse type is not 'Auto', and any parameter is non-zero then all other parameters must be non-zero. 
        """ 
        
        errors = []

        for pulse_key in ['Exc','Ref','Pump']:
            pulse_type = self.AdvPulseSelect[f'{pulse_key}Type'].currentText()
            if pulse_type != 'Auto':
                freq = self.AdvPulseSelect[f'{pulse_key}Freq'].value()
                length = self.AdvPulseSelect[f'{pulse_key}Length'].value()
                bandwidth = self.AdvPulseSelect[f'{pulse_key}Bandwidth'].value()
                params = [ length] # Frquency is allowed to be zero
                if pulse_type in ['Chirp','HS']:
                    params.append(bandwidth)
                non_zero_params = [p for p in params if np.abs(p) > 0]
                if 0 < len(non_zero_params) < len(params):
                    errors.append(f"For {pulse_key} pulse of type {pulse_type}, all parameters must be set if any are non-zero.")
        return errors


    def save_advaced_options(self):
        """
        Save the advanced options to a .yaml file
        """
        # Read current advanced options
        self.read_advanced_sequence_options()
        self.read_advanced_pulse_options()

        save_dict = {}
        save_dict['autoDEER_Vesion'] = ad.__version__
        save_dict['SequenceOptions'] = self.AdvSeqOptions
        save_dict['PulseOptions'] = self.AdvPulseOptions

        # Open a file dialog to save the options
        filename, _= QFileDialog.getSaveFileName(
            self,"Save Advanced Options", self.current_folder,"YAML (*.yaml)")
        
        if filename:
            path = Path(filename)
            filename_edit = str(path)
            with open(filename_edit, 'w') as f:
                yaml.dump(save_dict, f)

    def load_advanced_options(self):
        """
        Load advanced options from a .yaml file
        """
        filename, _= QFileDialog.getOpenFileName(
            self,"Load Advanced Options", self.current_folder,"YAML (*.yaml)")
        
        if filename:
            path = Path(filename)
            filename_edit = str(path)
            with open(filename_edit, 'r') as f:
                load_dict = yaml.safe_load(f)
        
        # Load Sequence Options
        if 'SequenceOptions' in load_dict:
            seq_options = load_dict['SequenceOptions']
            self.Exp_types.setCurrentText(seq_options['ExpType'].capitalize())
            self.update_advanced_sequence_options()
            for key, value in seq_options.items():
                if key == 'ExpType':
                    continue
                # Find the corresponding spinbox in SequenceParamRow
                for i in range(self.SequenceParamRow.count()):
                    item = self.SequenceParamRow.itemAt(i)
                    widget = item.widget()
                    if isinstance(widget, QtWidgets.QSpinBox):
                        label = self.SequenceParamRow.itemAt(i-1).widget()
                        if label.text() == key:
                            widget.setValue(value)
        
        # Load Pulse Options
        if 'PulseOptions' in load_dict:
            pulse_options = load_dict['PulseOptions']
            for key, value in pulse_options.items():
                if key in self.AdvPulseSelect:
                    widget = self.AdvPulseSelect[key]
                    if isinstance(widget, QtWidgets.QComboBox):
                        widget.setCurrentText(value)
                    elif isinstance(widget, QtWidgets.QDoubleSpinBox):
                        widget.setValue(value)

    def set_spectrometer_connected_light(self, state):
        if state == 0:
            light_pixmap = QtGui.QPixmap('icons:Red.png')
        elif state == 1:
            light_pixmap = QtGui.QPixmap('icons:Green.png')
        elif state == 2:
            light_pixmap = QtGui.QPixmap('icons:Yellow.png')
        
        light_pixmap = light_pixmap.scaledToHeight(30)
        self.Connected_Light.setPixmap(light_pixmap)

    def load_folder(self,*args, folder_path=None):
        if folder_path is None:
            folder_path = str(QFileDialog.getExistingDirectory(self, "Select Directory",str(Path.home())))
        self.pathLineEdit.setText(folder_path)
        self.current_folder = folder_path

        setup_logs(self.current_folder)
        global main_log
        main_log = logging.getLogger('autoDEER')
        main_log.info(f"Loading folder {self.current_folder}")

        # create folder for saving data if it doesn't exist
        if not os.path.exists(os.path.join(self.current_folder,'autoDEER_data')):
            os.makedirs(os.path.join(self.current_folder,'autoDEER_data'))
        self.data_folder = os.path.join(self.current_folder,'autoDEER_data')

    
    def load_epr_file(self, store_location):

        filename, _= QFileDialog.getOpenFileName(
            self,"Select a File", self.current_folder,"Data (*.DTA *.mat *.h5)")
        
        if filename:
                path = Path(filename)
                filename_edit = str(path)

        dataset = epr.eprload(filename_edit)
        self.current_data[store_location] = dataset

    def load_spectrometer_config(self, filename=None):
        
        if self.current_folder is None:
            QMessageBox.about(self,'ERORR!', 
                              'Please open a folder first.')
            return None

        if filename is None:
            filename, _= QFileDialog.getOpenFileName(
                self,"Select a spectrometer configuration file", str(Path.home()),"Data (*.yaml)")
        
        if filename:
            path = Path(filename)
            filename_edit = str(path)
        else:
            return None
        
        if self.spectromterInterface is not None:
            self.spectromterInterface = None
            if hasattr(self,'modeTuneButton'):
                self.modeTuneDialog.close()
                self.Resonator_layout.removeWidget(self.modeTuneButton)

        

        with open(filename_edit, mode='r') as f:
            config = yaml.safe_load(f)
            self.config = config
        
        self.logDialog = LogDialog(self)
        self.logDialog.hide()
        self.logbutton.clicked.connect(self.logDialog.show)

        main_log.info(f"Loading config file {filename_edit}")

        try:
            loglevels = config['autoDEER']['logging']
            change_log_level(loglevels['autoDEER'],loglevels['interface'])
        except KeyError:
            pass
        

        spectrometer = config['Spectrometer']

        if spectrometer['Manufacturer'] == 'Bruker':
            if spectrometer['AWG']:
                model = 'Bruker_AWG'
            elif spectrometer['MPFU']:
                model = 'Bruker_MPFU'
        elif spectrometer['Manufacturer'] == 'ETH':
            model = 'ETH_AWG'
        elif spectrometer['Manufacturer'] == 'Dummy':
            model = 'Dummy'
        
        try:
            if model == 'Dummy':
                from pyepr.hardware.dummy import dummyInterface
                self.spectromterInterface = dummyInterface(filename_edit)
                self.spectromterInterface.savefolder = self.current_folder
                self.Bruker=False
                self.pump_pulses = [epr.RectPulse,epr.ChirpPulse,epr.HSPulse]
            elif model == 'ETH_AWG':
                from pyepr.hardware.ETH_awg import ETH_awg_interface
                self.spectromterInterface = ETH_awg_interface(self.config)
                self.spectromterInterface.savefolder = self.data_folder
                self.Bruker=False
                self.modeTuneDialog = ModeTune(self.spectromterInterface, gyro=self.gyro, threadpool=self.threadpool, current_folder=self.current_folder)
                self.modeTuneDialog.dataUpdated.connect(self.update_resonator_info)
                self.modeTuneButton = QPushButton('Mode Tune')
                self.formLayout_2.addWidget(self.modeTuneButton)
                self.modeTuneButton.clicked.connect(self.modeTuneDialog.show)
                self.pump_pulses = [epr.RectPulse,epr.ChirpPulse]

            elif model == 'Bruker_MPFU':
                from pyepr.hardware.Bruker_MPFU import BrukerMPFU
                self.spectromterInterface = BrukerMPFU(filename_edit)
                self.spectromterInterface.savefolder = self.current_folder
                self.Bruker=True
                self.pump_pulses = [epr.RectPulse]
            elif model == 'Bruker_AWG':
                from pyepr.hardware.Bruker_AWG import BrukerAWG
                self.spectromterInterface = BrukerAWG(filename_edit)
                self.spectromterInterface.savefolder = self.current_folder
                self.Bruker=True
                self.pump_pulses = [epr.RectPulse,epr.ChirpPulse]
        except ImportError as e:
            QMessageBox.about(self,'ERORR!', 
                              'The spectrometer interface could not be loaded!\n'+
                              'Please check that the correct packages are installed!\n'+
                              'See the documentation for more information.')
            print(e)
            main_log.error('The spectrometer interface could not be loaded!')
            return None


        resonator_list = list(self.config['Resonators'].keys())

        # Find resonators
        self.resonatorComboBox.clear()
        self.resonatorComboBox.addItems(self.config['Resonators'].keys())
        self.resonatorComboBox.currentIndexChanged.connect(self.select_resonator)
        self.fcDoubleSpinBox.valueChanged.connect(self.change_LO)

        if len(resonator_list) == 0:
            QMessageBox.about(self,'ERORR!', 'No resonators found in config file!')
            main_log.error('No resonators found in config file!')
            return None
        # Set LO to resonator central frequency
        self.select_resonator()

        self.AWG = self.config['Spectrometer']['AWG']

        # Set the wavefrom precision to 1/Sampling Frequency
        if 'Waveform Precision' in self.config['Spectrometer']:
            waveform_precision = self.config['Spectrometer']['Waveform Precision']
        else:
            waveform_precision = 1/self.config['Spectrometer']['Bridge']['Sample Freq']
        epr.set_waveform_precision(waveform_precision)
        main_log.debug(f"Setting waveform precision to {epr.get_waveform_precision()}")

        # Get user preferences
        try:
            self.cores = int(self.config['autoDEER']['cores'])
            self.q_DEER.cores = self.cores
            self.longDEER.cores = self.cores
        except:
            self.q_DEER.cores = 1
            self.longDEER.cores = 1

        # Extract DeerLab fit parameters
        if ('autoDEER' in self.config):
            if ('DeerLab' in self.config['autoDEER']):
                self.DL_params = self.config['autoDEER']['DeerLab']
            else:
                self.DL_params = {}
            
            if ('Min_tp' in self.config['autoDEER']):
                self.Min_tp = self.config['autoDEER']['Min_tp']

        self.q_DEER.DL_params = self.DL_params
        self.longDEER.DL_params = self.DL_params
        

        try:
            self.waveform_precision = spectrometer['Bridge']['Pulse dt']
        except KeyError:
            self.waveform_precision = 1 


    def detect_recent_experiments(self,worker=None):
        # To cut down on runtime the current folder is checked for recent experiments

        all_files = os.listdir(self.current_folder)
        all_files.sort()
        EFDS_files = []
        Resonator_files = []
        SRT_scan_files = []
        CP_data = []
        Ref1D_data = []
        Tm_data = []
        Ref2D_data = []
        for file in all_files:
            if not file.endswith('h5'):
                continue
            if 'EDFS' in file:
                EFDS_files.append(file)
            elif 'ResPro' in file:
                Resonator_files.append(file)
            elif 'reptime' in file:
                SRT_scan_files.append(file)
            elif 'CarrPurcellSequence' in file:
                CP_data.append(file)
            elif 'RefocusedEcho1DSequence' in file:
                Ref1D_data.append(file)
            elif 'T2RelaxationSequence' in file:
                Tm_data.append(file)
            elif 'RefocusedEcho2DSequence' in file:
                Ref2D_data.append(file)

        # TODO: Add dialog warning box with files that are being loaded and ask if they should be. 
        # Create a dialog to ask the user which files to skip
        if any([EFDS_files, Resonator_files, SRT_scan_files, CP_data, Ref1D_data, Tm_data,Ref2D_data]):
            dialog = QDialog(self)
            dialog.setWindowTitle("Detected Previous Experiments")
            layout = QVBoxLayout()
            dialog.setLayout(layout)
            
            label = QLabel("<h2>Select experiments to skip:</h2>")
            label.setTextFormat(QtCore.Qt.TextFormat.RichText)
            layout.addWidget(label)
            
            checkboxes = {}
            
            if EFDS_files:
                cb = QPushButton(f"Field Sweep ({EFDS_files[-1]})")
                cb.setCheckable(True)
                cb.setChecked(True)
                checkboxes['run_fsweep'] = cb
                layout.addWidget(cb)
            
            if Resonator_files:
                cb = QPushButton(f"Resonator Profile ({Resonator_files[-1]})")
                cb.setCheckable(True)
                cb.setChecked(True)
                checkboxes['run_respro'] = cb
                layout.addWidget(cb)
            
            if SRT_scan_files:
                cb = QPushButton(f"Shot Repetition Time ({SRT_scan_files[-1]})")
                cb.setCheckable(True)
                cb.setChecked(True)
                checkboxes['run_reptime_opt'] = cb
                layout.addWidget(cb)
            
            if CP_data:
                cb = QPushButton(f"Carr-Purcell Relaxation ({CP_data[-1]})")
                cb.setCheckable(True)
                cb.setChecked(True)
                checkboxes['run_CP_relax'] = cb
                layout.addWidget(cb)
            
            if Ref1D_data:
                cb = QPushButton(f"Refocused Echo 1D ({Ref1D_data[-1]})")
                cb.setCheckable(True)
                cb.setChecked(True)
                checkboxes['run_1D_refocused_echo'] = cb
                layout.addWidget(cb)
            
            if Tm_data:
                cb = QPushButton(f"T2 Relaxation ({Tm_data[-1]})")

                cb.setCheckable(True)
                cb.setChecked(True)
                checkboxes['run_T2_relax'] = cb
                layout.addWidget(cb)
            
            if Ref2D_data:
                cb = QPushButton(f"Refocused Echo 2D ({Ref2D_data[-1]})")
                cb.setCheckable(True)
                cb.setChecked(True)
                checkboxes['run_2D_refocused_echo'] = cb
                layout.addWidget(cb)
            
            buttons = QHBoxLayout()
            ok_button = QPushButton("Skip Selected")
            ok_button.clicked.connect(dialog.accept)
            cancel_button = QPushButton("Skip None")
            cancel_button.clicked.connect(dialog.reject)
            buttons.addWidget(cancel_button)
            buttons.addWidget(ok_button)
            layout.addLayout(buttons)
            
            result = dialog.exec()

                        
            if result == 1:
                skip_list = [key for key, checkbox in checkboxes.items() if checkbox.isChecked()]

                progress = QProgressDialog("Processing Skipped Experiments...", None, 0, len(skip_list), self)
                progress.setWindowModality(QtCore.Qt.WindowModality.WindowModal)
                progress.setValue(0)                
                # Load the selected data files
                if 'run_fsweep' in skip_list and EFDS_files:
                    self.update_fieldsweep(epr.eprload(os.path.join(self.current_folder, EFDS_files[-1])), threaded=False, skip_recalc_d0=True)
                    progress.setValue(progress.value() + 1)
                if 'run_respro' in skip_list and Resonator_files:
                    self.update_respro(epr.eprload(os.path.join(self.current_folder, Resonator_files[-1])), threaded=False)
                    progress.setValue(progress.value() + 1)
                if 'run_reptime_opt' in skip_list and SRT_scan_files:
                    self.update_reptime(epr.eprload(os.path.join(self.current_folder, SRT_scan_files[-1])), threaded=False)
                    progress.setValue(progress.value() + 1)
                if 'run_CP_relax' in skip_list and CP_data:
                    self.update_relax(epr.eprload(os.path.join(self.current_folder, CP_data[-1])), threaded=False)
                    progress.setValue(progress.value() + 1)
                if 'run_1D_refocused_echo' in skip_list and Ref1D_data:
                    self.update_relax(epr.eprload(os.path.join(self.current_folder, Ref1D_data[-1])), threaded=False)
                    progress.setValue(progress.value() + 1)
                if 'run_T2_relax' in skip_list and Tm_data:
                    self.update_relax(epr.eprload(os.path.join(self.current_folder, Tm_data[-1])), threaded=False)
                    progress.setValue(progress.value() + 1)
                # if 'run_2D_refocused_echo' in skip_list and Ref2D_data:
                #     self.update_relax(epr.eprload(os.path.join(self.current_folder, Ref2D_data[-1])), threaded=False)
                #     progress.setValue(progress.value() + 1)

                progress.close()
                if worker is not None:
                    self.worker.update_skip_list(skip_list)



        # skip_list = []
        # if EFDS_files != []:
        #     self.update_fieldsweep(epr.eprload(os.path.join(self.current_folder,EFDS_files[-1])),threaded=False,skip_recalc_d0=True)
        #     skip_list.append('run_fsweep')
        # if Resonator_files != []:
        #     self.update_respro(epr.eprload(os.path.join(self.current_folder,Resonator_files[-1])),threaded=False)
        #     skip_list.append('run_respro')
        # if SRT_scan_files != []:
        #     self.update_reptime(epr.eprload(os.path.join(self.current_folder,SRT_scan_files[-1])),threaded=False)
        #     skip_list.append('run_reptime_opt')
        # if CP_data != []:
        #     self.update_relax(epr.eprload(os.path.join(self.current_folder,CP_data[-1])),threaded=False)
        #     skip_list.append('run_CP_relax')
        # if Ref1D_data != []:
        #     self.update_relax(epr.eprload(os.path.join(self.current_folder,Ref1D_data[-1])),threaded=False)
        #     skip_list.append('run_1D_refocused_echo')
        # if Tm_data != []:
        #     self.update_relax(epr.eprload(os.path.join(self.current_folder,Tm_data[-1])),threaded=False)
        #     skip_list.append('run_T2_relax')
        # if worker is not None:
        #     self.worker.update_skip_list(skip_list)

         


    def select_resonator(self):
        key = self.resonatorComboBox.currentText()
        main_log.info(f"Selecting resonator {key}")
        self.LO = self.config['Resonators'][key]['Center Freq']
        self.Q = self.config['Resonators'][key]['Q']
        self.fcDoubleSpinBox.setValue(self.LO)
        main_log.info(f"Setting LO to {self.LO} GHz")

    def update_resonator_info(self, kwargs):
        if 'fc' in kwargs:
            self.fcDoubleSpinBox.setValue(kwargs['fc'])
            self.LO = kwargs['fc']
            main_log.info(f"Setting LO to {self.LO} GHz")
        if 'q' in kwargs:
            self.qDoubleSpinBox.setValue(kwargs['q'])
            self.Q = kwargs['q']
            main_log.info(f"Setting Q to {kwargs['q']}")
    
    def change_LO(self):
        self.LO = self.fcDoubleSpinBox.value()
        main_log.info(f"Setting LO to {self.LO} GHz")

    @property
    def remaining_time(self):
        """
        Returns the remaining time in hours
        """
        max_time = self.userinput.get('MaxTime', 0)
        return max_time - (time.time() - self.starttime) / 3600

    def connect_spectrometer(self):
        if self.spectromterInterface is None:
            QMessageBox.about(self,'ERORR!', 'A interface needs to be loaded first!')
            main_log.error('Could not connect to spectrometer. A interface needs to be loaded first!')
            return None


        try:
            self.spectromterInterface.connect()
        except RuntimeError as e:
            QMessageBox.about(self, 'Connection ERORR!',str(e))
            main_log.error(f"Could not connect to spectrometer. {e}")
            return None
            

        self.connected = True
        self.set_spectrometer_connected_light(1)

    def raise_warning(self, message):
        # QDialog.about(self,'Warning!', message)
        main_log.warning(message)

    def save_data(self,dataset,experiment=None,folder='data'):
        if experiment is None:
            experiment = dataset.seq_name
        # add user_inputs to attributes
        # dataset.attrs.update(self.userinput)

        filename = create_save_name(self.userinput['sample'],experiment,True,self.userinput['project'],self.userinput['comment'])
        filename += ".h5"
        if folder == 'data':
            dataset.to_netcdf(os.path.join(self.data_folder,filename),engine='h5netcdf',invalid_netcdf=True)
        else:
            dataset.to_netcdf(os.path.join(self.current_folder,filename),engine='h5netcdf',invalid_netcdf=True)
        main_log.debug(f"Saving dataset to {filename}")



    def update_fieldsweep(self, dataset=None, threaded=True,**kwargs):

        if dataset is None:
            dataset = self.current_data['fieldsweep']
        else:
            self.current_data['fieldsweep'] = dataset
            self.save_data(dataset,'EDFS',folder='main')

        fsweep_analysis = epr.FieldSweepAnalysis(dataset)
        fsweep_analysis.calc_gyro()
        main_log.info(f"Calculated gyro {fsweep_analysis.gyro*1e3:.3f} MHz/T")
        if self.worker is not None:
            self.worker.set_noise_mode(np.min([fsweep_analysis.calc_noise_level(SNR_target=10),20]))
        
        self.setup_plot_ax.cla()
        fsweep_analysis.plot(axs=self.setup_plot_ax,fig=self.setup_figure_canvas.figure)
        self.setup_figure_canvas.draw()
        self.Tab_widget.setCurrentIndex(1)

        # For Bruker recalculate d0
        skip_recalc_d0 = kwargs.get('skip_recalc_d0',False)
        # if self.Bruker and not skip_recalc_d0:
        #     main_log.info("Calculating d0 from Hahn Echo")
        #     B = self.LO/fsweep_analysis.gyro
        #     self.spectromterInterface.calc_d0_from_Hahn_Echo(B=B, freq=self.LO)

        if self.worker is not None:
            self.worker.update_gyro(fsweep_analysis.gyro)
        try:
            fit = self.config['autoDEER']['FieldSweep']['Fit']
        except KeyError:
            fit = False
        
        if threaded:
            worker = Worker(fieldsweep_fit, fsweep_analysis,fit)
            worker.signals.result.connect(self.refresh_fieldsweep_after_fit)
            
            self.threadpool.start(worker)
        else:
            fit_result = fieldsweep_fit(fsweep_analysis,fit)
            self.refresh_fieldsweep_after_fit(fit_result)


 
    def refresh_fieldsweep_after_fit(self, fitresult=None):

        if (fitresult is None) and 'fieldsweep' in self.current_results:
            fitresult = self.current_results['fieldsweep']
        elif fitresult is None:
            return None
        else:
            self.current_results['fieldsweep'] = fitresult
        self.gyro = fitresult.gyro
        
        if self.waitCondition is not None: # Wake up the runner thread
            self.waitCondition.wakeAll()
        
        self.setup_plot_ax.cla()
        fitresult.plot(axs=self.setup_plot_ax,fig=self.setup_figure_canvas.figure)
        self.setup_figure_canvas.draw()

        # Add fit resultss

        self.gyroSpinBox.setValue(fitresult.gyro*1e3)
        
        self.Tab_widget.setCurrentIndex(1)

        if not ((self.pulses is None) or (self.pulses == {})):
            self.update_optimise_pulses_figure()



    def update_respro(self, dataset=None, threaded=True):

        if dataset is None:
            dataset = self.current_data['respro']
        else:
            self.current_data['respro'] = dataset
            self.save_data(dataset,'ResPro',folder='main')
        
        f_lims = (self.config['Spectrometer']['Bridge']['Min Freq'], self.config['Spectrometer']['Bridge']['Max Freq'])
        if threaded:

            # worker = Worker(respro_process, dataset, f_axis,self.current_results['fieldsweep'], cores=self.cores)
            worker = Worker(respro_process, dataset,f_lims, self.current_results['fieldsweep'], cores=self.cores)

            worker.signals.result.connect(self.refresh_respro)

            self.threadpool.start(worker)
        else:
            self.refresh_respro(respro_process(dataset,f_lims, self.current_results['fieldsweep'], cores=self.cores),threaded=False)


    def create_setup_figure(self):
        fig, axs  = plt.subplots(1,1,figsize=(12.5, 6.28))
        self.setup_figure_canvas = FigureCanvas(fig)
        Navbar = NavigationToolbar2QT(self.setup_figure_canvas, self)
        Navbar.setMaximumHeight(24)
        self.respro_v_left.addWidget(self.setup_figure_canvas)
        self.respro_v_left.addWidget(Navbar)
        self.setup_plot_ax = axs


    def refresh_respro(self, *args,threaded=True):
        if (len(args)>0)and (len(args[0]) == 2):
            LO = self.LO
            fitresult:epr.ResonatorProfileAnalysis = args[0][0]
            self.LO = args[0][1]
            self.LO = fitresult.results.fc

            LO_sweep_width = fitresult.freqs.max() - fitresult.freqs.min()
            LO_shift = self.LO - LO
            if np.abs(LO_shift) > LO_sweep_width:
                if LO_shift > 0:
                    self.LO + LO_sweep_width/2
                else:
                    self.LO - LO_sweep_width/2

        self.current_results['respro'] = fitresult
        self.spectromterInterface.resonator = fitresult


        if self.worker is not None:
            print(f"New center frequency: {self.LO:.2f} GHz")
            main_log.info(f"Setting center frequency to {self.LO:.6f} GHz")

            if LO_shift > 0.1:
                main_log.info(f"New center frequency is more than 100 MHz from previous value, repeating resonator profile")
                self.worker.update_freq(self.LO,repeat=True)
                optimise_pulses=False
                
            elif ('fieldsweep' in self.current_results) and (np.abs(self.current_results['fieldsweep'].freq - self.LO) > 0.025): # check if the last last field sweep was measured close the centre frequency
                main_log.info(f"Centre frequency is more than 25 MHz from last field sweep, repeating field sweep")
                self.worker.repeat_fieldsweep()
                optimise_pulses= True
            else:    
                self.worker.update_freq(self.LO,repeat=False)
                optimise_pulses=True

            

        self.Tab_widget.setCurrentIndex(1)    
        self.resonatorProfileFigure()
        # Add fit results

        self.centreFrequencyDoubleSpinBox.setValue(fitresult.results.fc)
        self.centreFrequencyCI.setText(f"({fitresult.results.fcUncert.ci(95)[0]:.2f},{fitresult.results.fcUncert.ci(95)[1]:.2f})")
        self.qDoubleSpinBox.setValue(fitresult.results.q)
        self.qCI.setText(f"({fitresult.results.qUncert.ci(95)[0]:.2f},{fitresult.results.qUncert.ci(95)[1]:.2f})")
        self.nu1doubleSpinBox.setValue(fitresult.model.max()*1e3)
        
        if optimise_pulses:
            main_log.info(f"Resonator centre frequency {fitresult.results.fc:.4f} GHz")
            self.optimise_pulses_button(threaded)
        else:
            if self.waitCondition is not None: # Wake up the runner thread
                self.waitCondition.wakeAll()


    def resonatorProfileFigure(self):
        if 'respro' in self.current_results:
            fitresult = self.current_results['respro']
        else:
            main_log.error('No resonator profile data found')
            return None

        self.setup_plot_ax.cla()

        if 'fieldsweep'in self.current_results:
            fitresult.plot(fieldsweep=self.current_results['fieldsweep'],axs=self.setup_plot_ax,fig=self.setup_figure_canvas.figure);
        else:
            fitresult.plot(axs=self.setup_plot_ax,fig=self.setup_figure_canvas.figure)

        self.setup_figure_canvas.draw()



    def optimise_pulses_button(self,threaded=True):
        if (self.pulses is None) or (self.pulses == {}):
            self.optimise_pulses(threaded=threaded)
        elif self.pulses['pump_pulse'] is None:
            self.optimise_pulses(threaded=threaded)
        else:
            self.update_optimise_pulses_figure()
        

    def optimise_pulses(self,threaded=True):
        if threaded:
            worker = Worker(self._optimise_pulses_in_background)
            worker.signals.result.connect(self.update_pulses)
            self.threadpool.start(worker)
        else:
            self.update_pulses(self._optimise_pulses_in_background())


    def _optimise_pulses_in_background(self):
        resonator = self.current_results['respro']
        spectrum = self.current_results['fieldsweep']
        AdvPulse_types = {}
        skip_optimisation = False

        if self.AdvPulseOptions != {}: # If advance pulse options exist
            # If all pulse types are auto, then continue as normal
            all_auto = True
            all_zero = True
            
            for pulse_key in ['Exc','Ref','Pump']:
                pulse_type = self.AdvPulseOptions.get(f'{pulse_key}Type','Auto')
                # Remove any "." from pulse type
                pulse_type = pulse_type.replace('.','')
                pulse_tp = self.AdvPulseOptions.get(f'{pulse_key}Length',0)
                if pulse_type.lower() != 'auto':
                    all_auto = False
                    if pulse_type.lower() == 'rect':
                        AdvPulse_types[f"{pulse_key}PulseShape"] = epr.RectPulse
                    elif pulse_type.lower() == 'chirp':
                        AdvPulse_types[f"{pulse_key}PulseShape"] = epr.ChirpPulse
                    elif pulse_type.lower() == 'hs':
                        AdvPulse_types[f"{pulse_key}PulseShape"] = epr.HSPulse
                    elif pulse_type.lower() == 'gauss':
                        AdvPulse_types[f"{pulse_key}PulseShape"] = epr.GaussianPulse
                    else:
                        raise ValueError(f"Unknown pulse type {pulse_type} for {pulse_key} pulse")
                if pulse_tp > 0:
                    all_zero = False

            # Create all pulses:
            if not all_auto and not all_zero:
                pulses = {}
                skip_optimisation = True
                main_log.info(f"Creating pulses from user defined parameters")
                for pulse_key in ['Exc','Ref','Pump']:
                    pulse_params = {}
                    pulse_type = AdvPulse_types[f"{pulse_key}PulseShape"]
                    pulse_params['tp'] = self.AdvPulseOptions.get(f'{pulse_key}Length')
                    pulse_params['freq'] = self.AdvPulseOptions.get(f'{pulse_key}Freq')
                    if issubclass(pulse_type,epr.FrequencySweptPulse):
                        pulse_params['BW'] = self.AdvPulseOptions.get(f'{pulse_key}Bandwidth')
                    if pulse_key == 'Exc':
                        pulse_params['flipangle'] = np.pi/2
                    else:
                        pulse_params['flipangle'] = np.pi
                    pulses[pulse_key.lower()+'_pulse'] = pulse_type(**pulse_params)
                pulses['det_event'] = epr.Detection(freq=pulses['exc_pulse'].freq, tp=2*pulses['exc_pulse'].tp.value)


        if "PumpPulseShape" in AdvPulse_types:
            self.pump_pulses = [AdvPulse_types.pop("PumpPulseShape")]
        if skip_optimisation:
            main_log.info(f"Creating pulses with user defined pulse types and lengths")
            self.userinput['AdvUserPulses'] = True
            self.userinput['AdvPulseShapes'] = False
        elif (not skip_optimisation) and (self.pulses == {}):  # No pulses have been created yet
            # self.pulses = ad.build_default_pulses(self.AWG,tp = self.Min_tp)
            if AdvPulse_types != {}:
                main_log.info(f"Creating pulses with user defined pulse types")
                self.userinput['AdvUserPulses'] = False
                self.userinput['AdvPulseShapes'] = False
            else:
                main_log.info(f"Creating pulses with auto defined pulse types")
                self.userinput['AdvUserPulses'] = False
                self.userinput['AdvPulseShapes'] = True  
            pulses = ad.create_pulses_shape(resonatorProfile=resonator,spectrum=spectrum, test_pulse_shapes=self.pump_pulses,**AdvPulse_types)
        
        elif (not skip_optimisation):  # Reoptimise pulses
            self.userinput['AdvUserPulses'] = False
            self.userinput['AdvPulseShapes'] = False
            if 'quickdeer' in self.current_results:
                ROI = self.current_results['quickdeer'].ROI
                r_min = ROI[0]
                ROI_width = ROI[1] - ROI[0]
                if (not ad.check_pulses_max_length(self.pulses.values(), r_min)) and (ROI_width < 4):
                    pulses = ad.create_pulses_shape(resonatorProfile=resonator,spectrum=spectrum,r_min=r_min, test_pulse_shapes=self.pump_pulses,**AdvPulse_types)
                    main_log.info(f"Pulse lenths adjusted, setting new pulses")
                elif (not ad.check_pulses_max_length(self.pulses.values(), r_min)):
                    main_log.info(f"ROI is detected to be broad, not adjusting pulse lengths. Potentially not optimal")
                    pulses  =self.pulses
                    return None
                else:
                    main_log.info(f"Pulse lenths are ok for ROI, not adjusting pulse lengths")
                    pulses  =self.pulses
                    return None
        
        est_lambda = ad.calc_est_modulation_depth(spectrum, **pulses, respro=resonator)
        main_log.info(f"Estimated modulation depth from pulses: {est_lambda:.3f}")
        return pulses, est_lambda
    


    def update_pulses(self,*args):
        if (args[0] is not None) and len(args[0]) == 2:
            pulses, est_lambda = args[0]
            self.pulses = pulses
            self.est_lambda = est_lambda
        
            pump_pulse = self.pulses['pump_pulse']
            ref_pulse = self.pulses['ref_pulse']
            exc_pulse = self.pulses['exc_pulse']

            
            # self.pulses = {'pump_pulse':pump_pulse, 'exc_pulse':exc_pulse, 'ref_pulse':ref_pulse}    
            self.pulses['pump_pulse'] = pump_pulse
            self.pulses['ref_pulse'] = ref_pulse
            self.pulses['exc_pulse'] = exc_pulse
            self.pulses['det_event'].freq = exc_pulse.freq
            self.pulses['det_event'].tp.value = 2*np.max([exc_pulse.tp.value,ref_pulse.tp.value])
            main_log.info(f"Optimised pulses")
            if self.worker is not None:
                self.worker.new_pulses(self.pulses)
            else:
                main_log.error(f"Can't send new pulses to worker as worker not found ")
        else:
            main_log.error("Incorrect inputs to update pulses")
        
        if self.waitCondition is not None: # Wake up the runner thread
            self.waitCondition.wakeAll()

        self.update_optimise_pulses_figure()
        self.Tab_widget.setCurrentIndex(1)

    
    def update_optimise_pulses_figure(self):
        try:
            pump_pulse = self.pulses['pump_pulse']
            exc_pulse = self.pulses['exc_pulse']
            ref_pulse = self.pulses['ref_pulse']
        except KeyError:
            main_log.error('No pulses found')
            return None

        self.setup_plot_ax.cla()
        ad.plot_overlap(self.current_results['fieldsweep'], pump_pulse, exc_pulse,ref_pulse, axs=self.setup_plot_ax,fig=self.setup_figure_canvas.figure, respro=self.current_results['respro'])
        self.setup_figure_canvas.draw()

        # update the pulse parameter grid
        type_to_pulse_hash = {epr.RectPulse:'Rect', epr.ChirpPulse:'Chirp', epr.HSPulse:'HS'}
        if isinstance(exc_pulse, epr.RectPulse):
            self.ExcFreqBox.setValue(param_in_MHz(exc_pulse.freq))
            est_BW = 1/(2*exc_pulse.tp.value) *1e3
            self.ExcBWBox.setValue(est_BW)
            # self.ExcBWBox.setSuffix(' ns')
            self.ExcLengthBox.setValue(exc_pulse.tp.value)
            self.ExcTypeLine.setText('Rect')
        else:
            center_freq = (param_in_MHz(exc_pulse.final_freq) + param_in_MHz(exc_pulse.init_freq))/2
            self.ExcFreqBox.setValue(center_freq)
            self.ExcBWBox.setValue(param_in_MHz(exc_pulse.bandwidth))
            # self.ExcBWBox.setSuffix(' MHz')
            self.ExcLengthBox.setValue(exc_pulse.tp.value)
            self.ExcTypeLine.setText(type_to_pulse_hash[type(exc_pulse)])
        
        if isinstance(ref_pulse, epr.RectPulse):
            self.RefFreqBox.setValue(param_in_MHz(ref_pulse.freq))
            est_BW = 1/(2*ref_pulse.tp.value) *1e3
            self.RefBWBox.setValue(est_BW)
            # self.RefBWBox.setSuffix(' ns')
            self.RefLengthBox.setValue(ref_pulse.tp.value)
            self.RefTypeLine.setText('Rect')
        else:
            center_freq = (param_in_MHz(ref_pulse.final_freq) + param_in_MHz(ref_pulse.init_freq))/2
            self.RefFreqBox.setValue(center_freq)
            self.RefBWBox.setValue(param_in_MHz(ref_pulse.bandwidth))
            # self.RefBWBox.setSuffix(' MHz')
            self.RefLengthBox.setValue(ref_pulse.tp.value)
            self.RefTypeLine.setText(type_to_pulse_hash[type(ref_pulse)])
        
        if isinstance(pump_pulse, epr.RectPulse):
            self.PumpFreqBox.setValue(param_in_MHz(pump_pulse.freq))
            est_BW = 1/(2*pump_pulse.tp.value) *1e3
            self.PumpBWBox.setValue(est_BW)
            # self.PumpBWBox.setSuffix(' ns')
            self.PumpLengthBox.setValue(pump_pulse.tp.value)
            self.PumpTypeLine.setText('Rect')
        else:
            center_freq = (param_in_MHz(pump_pulse.final_freq) + param_in_MHz(pump_pulse.init_freq))/2
            self.PumpFreqBox.setValue(center_freq)
            self.PumpBWBox.setValue(param_in_MHz(pump_pulse.bandwidth))
            # self.PumpBWBox.setSuffix(' MHz')
            self.PumpLengthBox.setValue(pump_pulse.tp.value)
            self.PumpTypeLine.setText(type_to_pulse_hash[type(pump_pulse)])

    def create_relax_figure(self,n_plots=3):
        # if hasattr(self,'relax_canvas'):
        #     self.relax_v_left :QVBoxLayout
        #     for i in reversed(range(self.relax_v_left.count())): 
        #         self.relax_v_left.itemAt(i).widget().setParent(None)

        if not hasattr(self,'relax_canvas'):
            fig  = plt.figure(figsize=(8,8), layout='constrained')
            self.relax_canvas = FigureCanvas(fig)
            Navbar = NavigationToolbar2QT(self.relax_canvas, self)
            Navbar.setMaximumHeight(24)
            self.relax_v_left.addWidget(self.relax_canvas)
            self.relax_v_left.addWidget(Navbar)
        else:
            fig = self.relax_canvas.figure

        fig.clear()

        gs = GridSpec(2, 2, figure=fig)
        if n_plots == 3:
            ax1 = fig.add_subplot(gs[0, :])
            ax2 = fig.add_subplot(gs[1, 0])
            ax3 = fig.add_subplot(gs[1, 1])

            self.relax_ax = [ax1, ax2, ax3]
        elif n_plots == 4:
            ax1 = fig.add_subplot(gs[0, 0]) # 1D relax
            ax2 = fig.add_subplot(gs[1, 0]) # 2D relax
            ax3 = fig.add_subplot(gs[0, 1]) # Reptime
            ax4 = fig.add_subplot(gs[1, 1]) # Optimal Tau

            self.relax_ax = [ax1, ax2, ax3, ax4]

    def refresh_relax_figure(self):
    
        fig = self.relax_canvas.figure
        self.relax_ax[0].cla()
        relax1D_results = []
        if 'CP-relax' in self.current_results:
            relax1D_results.append(self.current_results['CP-relax'])
        if 'Tm-relax' in self.current_results:
            relax1D_results.append(self.current_results['Tm-relax'])
        if 'RefEcho1D' in self.current_results:
            relax1D_results.append(self.current_results['RefEcho1D'])

        epr.plot_1Drelax(*relax1D_results, axs=self.relax_ax[0], fig=fig, cmap=ad.primary_colors)
            
        if 'RefEcho2D' in self.current_results:
            self.relax_ax[1].cla()
            self.current_results['RefEcho2D'].plot2D(axs=self.relax_ax[1], fig=fig)

        self.relax_canvas.draw()

    def initialise_deer_settings(self):
        

        if self.est_lambda is None:
            self.est_lambda = 0.4 
        self.aim_time = 2
        self.aim_MNR = 20

        # Assemble all relaxation data
        if (self.Exp_types.currentText() == '4pDEER'):
            exp = '4pDEER'
        else:
            exp = 'auto'

        if self.userinput['priority'].lower() == 'single':
            exp = '4pDEER'
            aim_SNR = self.priorties[self.userinput['priority']]
            self.aim_time = self.MaxTime.value() - ((time.time() - self.starttime) / (60*60)) # in hours
            self.aim_MNR = aim_SNR
        else:
            aim_SNR = self.aim_MNR/(self.est_lambda*self.label_eff)
            

        relax_data = {}
        if 'CP-relax' in self.current_results:
            relax_data['CP'] = self.current_results['CP-relax']
        if 'RefEcho1D' in self.current_results:
            relax_data['RefEcho1D'] = self.current_results['RefEcho1D']
        if 'Tm-relax' in self.current_results:
            relax_data['Tm'] = self.current_results['Tm-relax']
        if 'RefEcho2D' in self.current_results:
            relax_data['RefEcho2D'] = self.current_results['RefEcho2D']    
        
        #debug only, remove later
        store_pickle(relax_data,os.path.join(self.data_folder,'relax_data.pkl'))
        main_log.info(f'Calculating DEER settings with aim SNR {aim_SNR:.2f}, aim time {self.aim_time}hrs')
        self.deer_settings = ad.calc_DEER_settings(relax_data,exp, self.aim_time, aim_SNR,self.waveform_precision)
                
        main_log.debug('Calculated DEER Settings')
        self.deer_settings['criteria'] = self.aim_MNR

        self.deer_settings['autoStop'] = self.Time_autoStop_checkbox.isChecked()
        main_log.debug('Updating DEER Settings')
        
        # self.worker.update_deersettings(self.deer_settings)
        self.worker.signals.update_deer_settings.emit(self.deer_settings)
        main_log.debug('Updating Figure')
        self.update_tau_delays_figure([aim_SNR],[self.aim_time],labels=[f"MNR = {self.aim_MNR}"])
        self.OptimalExperiment.setText(self.deer_settings['ExpType'])
        if self.deer_settings['ExpType'] == '4pDEER':
            tau_evo = self.deer_settings['tau2']
        elif self.deer_settings['ExpType'] == '5pDEER':
            tau_evo = self.deer_settings['tau1'] + self.deer_settings['tau2']
        
        self.DipolarEvo2hrs.setValue(tau_evo)

        
        main_log.info(f"tau1 set to {self.deer_settings['tau1']:.2f} us")
        main_log.info(f"tau2 set to {self.deer_settings['tau2']:.2f} us")
        main_log.info(f"DEER Sequence set to {self.deer_settings['ExpType']}")

        return self.deer_settings
    
    def update_deer_settings(self, update_pulses=True, remaining_time=None, threaded=False):
        """
        Update DEER settings based on the current results and user input.
        If remaining_time is provided, it will use that instead of the calculated remaining time.
        """
        if self.userinput['priority'].lower() == 'single':
            func = self._update_single_deer_settings_and_pulses
        else:
            func = self._update_deer_settings_and_pulses
        
        main_log.info("Updating DEER settings")
        if threaded:
            if not update_pulses:
                raise ValueError("update_pulses must be True when threaded is True")
            worker = Worker(func, remaining_time)
            worker.signals.result.connect(self.update_pulses)
            self.threadpool.start(worker)
        else:
            if update_pulses:
                self.update_pulses(func(remaining_time,update_pulses))
            else:
                func(remaining_time,update_pulses)
        

    def _update_single_deer_settings_and_pulses(self, remaining_time=None,update_pulses=False):
        """
        Update DEER settings for single mode based on the current results and user input.
        If remaining_time is provided, it will use that instead of the calculated remaining time.
        """
        main_log.info("Updating DEER settings for single mode")

        if self.userinput['priority'].lower() == 'single':
            SNR_target = self.priorties[self.userinput['priority']]
            MNR_target = SNR_target
            single_mode = True
            exp = '4pDEER'

        relax_data = {}
        if 'CP-relax' in self.current_results:
            relax_data['CP'] = self.current_results['CP-relax']
        if 'RefEcho1D' in self.current_results:
            relax_data['Tm'] = self.current_results['RefEcho1D']
        elif 'Tm-relax' in self.current_results:
            relax_data['Tm'] = self.current_results['Tm-relax']
        if 'RefEcho2D' in self.current_results:
            relax_data['Ref2D'] = self.current_results['RefEcho2D']
                 
        #debug only, remove later
        store_pickle(relax_data,os.path.join(self.data_folder,'relax_data.pkl'))

        self.deer_settings = ad.calc_DEER_settings(relax_data,exp,remaining_time,SNR_target,self.waveform_precision,corr_factor=self.correction_factor)
        
        self.deer_settings['criteria'] = MNR_target
        
        self.deer_settings['autoStop'] = self.Time_autoStop_checkbox.isChecked()
        # self.worker.update_deersettings(self.deer_settings)
        self.worker.signals.update_deer_settings.emit(self.deer_settings)
        self.update_tau_delays_figure([SNR_target],[remaining_time],labels=[f"MNR = {MNR_target}"])
        
        main_log.info(f"tau1 set to {self.deer_settings['tau1']:.2f} us")
        main_log.info(f"tau2 set to {self.deer_settings['tau2']:.2f} us")
        main_log.info(f"DEER Sequence set to {self.deer_settings['ExpType']}")
        

    def _update_deer_settings_and_pulses(self,remaining_time=None,update_pulses=True):

        relax_data = {}
        if 'CP-relax' in self.current_results:
            relax_data['CP'] = self.current_results['CP-relax']
        if 'RefEcho1D' in self.current_results:
            relax_data['Tm'] = self.current_results['RefEcho1D']
        elif 'Tm-relax' in self.current_results:
            relax_data['Tm'] = self.current_results['Tm-relax']
        if 'RefEcho2D' in self.current_results:
            relax_data['Ref2D'] = self.current_results['RefEcho2D']
        
        if 'quickdeer' in self.current_results:
            data = self.current_results['quickdeer']
            rec_tau = self.current_results['quickdeer'].rec_tau_max
            dt = self.current_results['quickdeer'].rec_dt * 1e3
            dt = epr.round_step(dt,self.waveform_precision)
            dt= 8
            mod_depth = data.MNR * data.noiselvl
            if mod_depth > 1.0:
                mod_depth = 1.0

            mod_depth_correction = mod_depth/ self.est_lambda
            main_log.info(f"Modulation depth correction factor {mod_depth_correction:.3f}")

            if self.deer_settings['ExpType'] == '4pDEER':
                relax = self.current_results['Tm-relax']
            elif self.deer_settings['ExpType'] == '5pDEER':
                relax = self.current_results['CP-relax']

            self.correction_factor = ad.calc_correction_factor(relax,self.current_results['quickdeer'])

            if remaining_time is None:
                remaining_time = self.remaining_time
                main_log.debug(f"Remaining time {remaining_time:.2f} hours")
            else:
                remaining_time = remaining_time
                main_log.debug(f"Measuring DEER for {remaining_time:.2f} hours")

            MNR_target = self.priorties[self.userinput['priority']]

            ROI = self.current_results['quickdeer'].ROI
            r_min = ROI[0]
            ROI_width = ROI[1] - ROI[0]

            if (ROI_width < 4) and (r_min < 3.5):
                r_min = ROI[0]
                main_log.info(f"Detected narrow distance distribution, setting r_min to {r_min} nm")
            else: 
                r_min = 3.5
                main_log.info(f"Detected broad distance distribution, setting r_min to default {r_min} nm")

        else:
            dt = 8
            if self.est_lambda is None:
                self.est_lambda = 0.4 
            self.aim_time = 2
            self.aim_MNR = 20
            MNR_target = self.aim_MNR
            mod_depth_correction = self.label_eff
            remaining_time = self.aim_time

            r_min = 3.5
            rec_tau = None

        # Assemble all relaxation data
        if (self.Exp_types.currentText() == '4pDEER'):
            exp = '4pDEER'
        else:
            exp = 'auto'

        EDFS_analysis = self.current_results['fieldsweep']
        ResProAnalysis = self.current_results['respro']

        # Don't update pulses if they were manually set
        if self.userinput.get('AdvUserPulses',False):
            update_pulses = False
        
        # Check if the any pulses shapes were manually set
        AdvPulse_types = {}
        if self.userinput.get('AdvPulseShapes',False):
            for pulse_key in ['Exc','Ref','Pump']:
                pulse_type = self.AdvPulseOptions.get(f'{pulse_key}Type','auto')
                pulse_type = pulse_type.replace('.','')
                if pulse_type.lower() != 'auto':
                    if pulse_type.lower() == 'rect':
                        AdvPulse_types[f"{pulse_key}PulseShape"] = epr.RectPulse
                    elif pulse_type.lower() == 'chirp':
                        AdvPulse_types[f"{pulse_key}PulseShape"] = epr.ChirpPulse
                    elif pulse_type.lower() == 'hs':
                        AdvPulse_types[f"{pulse_key}PulseShape"] = epr.HSPulse
                    elif pulse_type.lower() == 'gauss':
                        AdvPulse_types[f"{pulse_key}PulseShape"] = epr.GaussianPulse
                    else:
                        raise ValueError(f"Unknown pulse type {pulse_type} for {pulse_key} pulse")

         
        #debug only, remove later
        store_pickle(relax_data,os.path.join(self.data_folder,'relax_data.pkl'))

        if self.fixed_tau is not None:
            # Using advanced mode with fixed tau
            if remaining_time is None:
                remaining_time = self.remaining_time
                main_log.debug(f"Remaining time {remaining_time:.2f} hours")
            else:
                remaining_time = remaining_time
                main_log.debug(f"Measuring DEER for {remaining_time:.2f} hours")
            exp = 'adv'
            MNR_target = self.priorties[self.userinput['priority']]
            mod_depth_correction = self.label_eff
            r_min = 3.5
            rec_tau = None

            if self.Exp_types.currentText() == 'auto' or self.Exp_types.currentText() == '5pDEER':
                main_log.info("Advanced mode selected with fixed tau, overriding auto/5pDEER/4pDEER selection")
                if update_pulses:
                    optimal_pulses_5p = ad.create_pulses_shape(ResProAnalysis,EDFS_analysis,n_pump_pulses=2,test_pulse_shapes=self.pump_pulses,verbosity=0,r_min=r_min,**AdvPulse_types)
                else:
                    optimal_pulses_5p = self.pulses
                functional_5p = ad.calc_functional(EDFS_analysis,**optimal_pulses_5p,resonator=ResProAnalysis,n_pump_pulses=2)
                mod_depth_5p = ad.calc_est_modulation_depth(EDFS_analysis,**optimal_pulses_5p,resonator=ResProAnalysis,n_pump_pulses=2)
                SNR_target = MNR_target/(mod_depth_5p * mod_depth_correction)
                deer_settings_5p = self.deer_settings
                deer_settings_4p = {}
                tau2_4p=0
                tau1_5p = deer_settings_5p['tau1']
                functional_4p = 0
                
                main_log.debug(f"5p DEER:  F= {functional_5p:.3f} tau_evo={deer_settings_5p['tau1'] + deer_settings_5p['tau2']:.3f}")
            elif self.Exp_types.currentText() == '4pDEER':
                main_log.info("Advanced mode selected with fixed tau, overriding auto/5pDEER/4pDEER selection")
                if update_pulses:
                    optimal_pulses_4p = ad.create_pulses_shape(ResProAnalysis,EDFS_analysis,n_pump_pulses=1,test_pulse_shapes=self.pump_pulses,verbosity=0,r_min=r_min,**AdvPulse_types) 
                else:
                    optimal_pulses_5p = self.pulses
                functional_4p = ad.calc_functional(EDFS_analysis,**optimal_pulses_4p,resonator=ResProAnalysis,n_pump_pulses=1)
                mod_depth_4p = ad.calc_est_modulation_depth(EDFS_analysis,**optimal_pulses_4p,resonator=ResProAnalysis,n_pump_pulses=1)
                SNR_target = MNR_target/(mod_depth_4p * mod_depth_correction)
                deer_settings_4p = self.deer_settings
                deer_settings_5p = {}
                tau1_5p = 0
                tau2_4p = deer_settings_4p['tau2']
                functional_5p = 0
                main_log.debug(f"4p DEER:  F= {functional_4p:.3f} tau_evo={tau2_4p:.3f}")
        
        else:

            # Check 4p DEER with 1 pump vs 5p DEER with 2 pumps
            if exp == '4pDEER' or exp =='auto':
                if update_pulses:
                    optimal_pulses_4p = ad.create_pulses_shape(ResProAnalysis,EDFS_analysis,n_pump_pulses=1,test_pulse_shapes=self.pump_pulses,verbosity=0,r_min=r_min,**AdvPulse_types)
                else:
                    optimal_pulses_4p = self.pulses
                functional_4p = ad.calc_functional(EDFS_analysis,**optimal_pulses_4p,resonator=ResProAnalysis,n_pump_pulses=1)
                mod_depth_4p = ad.calc_est_modulation_depth(EDFS_analysis,**optimal_pulses_4p,resonator=ResProAnalysis,n_pump_pulses=1)

                SNR_target = MNR_target/(mod_depth_4p * mod_depth_correction)
                main_log.info(f"4p DEER, mod_depth = {mod_depth_4p:.2f}, SNR target {SNR_target:.2f}")
                deer_settings_4p = ad.calc_DEER_settings(relax_data,'4pDEER',remaining_time,SNR_target,self.waveform_precision,corr_factor=self.correction_factor,rec_tau=rec_tau)
                tau2_4p = deer_settings_4p['tau2']
                main_log.debug(f"4p DEER:  F= {functional_4p:.3f} tau_evo={tau2_4p:.3f}")

            else:
                deer_settings_4p = {}
                tau2_4p=0
                functional_4p = 0

            if exp == '5pDEER' or exp == 'auto':
                if update_pulses:
                    optimal_pulses_5p = ad.create_pulses_shape(ResProAnalysis,EDFS_analysis,n_pump_pulses=2,test_pulse_shapes=self.pump_pulses,verbosity=0,r_min=r_min,**AdvPulse_types)
                else:
                    optimal_pulses_5p = self.pulses
                functional_5p = ad.calc_functional(EDFS_analysis,**optimal_pulses_5p,resonator=ResProAnalysis,n_pump_pulses=2)
                mod_depth_5p = ad.calc_est_modulation_depth(EDFS_analysis,**optimal_pulses_5p,resonator=ResProAnalysis,n_pump_pulses=2)

                SNR_target = MNR_target/(mod_depth_5p * mod_depth_correction)
                main_log.info(f"5p DEER, mod_depth = {mod_depth_5p:.2f}, SNR target {SNR_target:.2f}")
                deer_settings_5p = ad.calc_DEER_settings(relax_data,'5pDEER',remaining_time,SNR_target,self.waveform_precision,corr_factor=self.correction_factor,rec_tau=rec_tau)
                tau1_5p = deer_settings_5p['tau1']
                main_log.debug(f"5p DEER:  F= {functional_5p:.3f} tau_evo={deer_settings_5p['tau1'] + deer_settings_5p['tau2']:.3f}")

            else:
                deer_settings_5p = {}
                tau1_5p = 0
                functional_5p = 0
        
        if deer_settings_4p == {} and deer_settings_5p == {}:
            main_log.error("No DEER settings found")
            return None

        # Select between five-pulse and four-pulse DEER
        if (tau2_4p >= 2*tau1_5p) and (functional_4p > functional_5p * 0.9):
            autoDEERmode = '4pDEER'
        else:
            autoDEERmode = '5pDEER'
        if (tau1_5p <1) and (tau2_4p > 0.1):
            autoDEERmode = '4pDEER'
        main_log.info(f"AutoDEER mode set to {autoDEERmode}")

        if autoDEERmode == '4pDEER':
            self.deer_settings = deer_settings_4p
            new_pulses = optimal_pulses_4p
            new_est_lambda = mod_depth_4p * mod_depth_correction
        elif autoDEERmode == '5pDEER':
            self.deer_settings = deer_settings_5p
            new_pulses = optimal_pulses_5p
            new_est_lambda = mod_depth_5p * mod_depth_correction
            

        # self.deer_settings = ad.calc_DEER_settings(relax_data,exp,remaining_time,SNR_target,self.waveform_precision,corr_factor=self.correction_factor,rec_tau=rec_tau)

        # self.deer_settings['dt'] = dt
        self.deer_settings['criteria'] = MNR_target
        
        self.deer_settings['autoStop'] = self.Time_autoStop_checkbox.isChecked()
        # self.worker.update_deersettings(self.deer_settings)
        self.worker.signals.update_deer_settings.emit(self.deer_settings)
        self.update_tau_delays_figure([SNR_target],[remaining_time],labels=[f"MNR = {MNR_target}"])
        
        main_log.info(f"tau1 set to {self.deer_settings['tau1']:.2f} us")
        main_log.info(f"tau2 set to {self.deer_settings['tau2']:.2f} us")
        main_log.info(f"DEER Sequence set to {self.deer_settings['ExpType']}")

        if update_pulses:
            return new_pulses, new_est_lambda
        else:
            self.est_lambda = new_est_lambda
            return self.deer_settings

    def update_tau_delays_figure(self, SNRs, MeasTimes, labels=None):

        fig = self.relax_canvas.figure
        axs = self.relax_ax[-1]
        axs.cla()
        
        if 'CP-relax' in self.current_results:
            CP_analysis = self.current_results['CP-relax']
            ad.plot_optimal_tau(CP_analysis,SNRs,MeasTimes,MaxMeasTime=36, labels=['5pDEER'],fig=fig,axs=axs,cmap=[epr.primary_colors[0]],corr_factor=self.correction_factor)

        if 'RefEcho2D' in self.current_results:
            main_log.debug('Using RefEcho2D for optimal tau calculation')
            ad.plot_optimal_tau(self.current_results['RefEcho2D'],SNRs,MeasTimes,MaxMeasTime=36, labels=['4pDEER'],fig=fig,axs=axs,cmap=[epr.primary_colors[1]],corr_factor=self.correction_factor)
        elif 'RefEcho1D' in self.current_results:
            main_log.debug('Using RefEcho1D for optimal tau calculation')
            ad.plot_optimal_tau(self.current_results['RefEcho1D'],SNRs,MeasTimes,MaxMeasTime=36, labels=['4pDEER'],fig=fig,axs=axs,cmap=[epr.primary_colors[1]],corr_factor=self.correction_factor)
        elif 'Tm-relax' in self.current_results:
            main_log.debug('Using Tm-relax for optimal tau calculation')
            ad.plot_optimal_tau(self.current_results['Tm-relax'],SNRs,MeasTimes,MaxMeasTime=36, labels=['4pDEER'],fig=fig,axs=axs,cmap=[epr.primary_colors[1]],corr_factor=self.correction_factor)

        axs.set_title(labels[0])

    def update_relax2D(self, dataset=None):
        if dataset is None:
            dataset = self.current_data['RefEcho2D']
        else:
            self.current_data['RefEcho2D'] = dataset
            self.save_data(dataset,folder='main')

        # Since there is no fitting in the 2D data analysis it can be run in the main threads

        relax2DAnalysis = ad.RefocusedEcho2DAnalysis(dataset)
        self.current_results['RefEcho2D'] = relax2DAnalysis

        self.refresh_relax_figure()

        self.update_deer_settings(update_pulses=False, threaded=False)

        if self.waitCondition is not None: # Wake up the runner thread
            self.waitCondition.wakeAll()
        
    def update_relax(self, dataset, threaded=True):
        # Recieves dataset from the worker thread and updates the relax figure
        
        dataset = dataset.epr.correctphasefull # Correct the phase of the dataset

        short_name = short_name_dict[dataset.seq_name]

        
        # check if the relaxation data is being merged

        if short_name in self.current_data:
            # attempt to merge datasets only if the new first datapoint is > than the last old datapoint else replace the old dataset
            if 'tau' in dataset.coords:
                axis_label = 'tau'
            elif 't' in dataset.coords:
                axis_label = 't'
            elif 'tau_1' in dataset.coords:
                axis_label = 'tau_1'
            elif 'tau_2' in dataset.coords:
                axis_label = 'tau_2'
            elif "tau1" in dataset.coords:
                axis_label = 'tau1'
            elif "tau2" in dataset.coords:
                axis_label = 'tau2'
            else:
                axis_label = 'X'

            new_dataset_first = getattr(dataset,axis_label).values[0]
            old_dataset_last = getattr(self.current_data[short_name],axis_label).values[-1]
            main_log.debug(f"New dataset first {new_dataset_first:.2f}, old dataset last {old_dataset_last:.2f}")

            if new_dataset_first >= (old_dataset_last - 1e-2):  # Allows for a small error
                main_log.info(f'Merging {short_name} datasets')
                dataset = self.current_data[short_name].epr.merge(dataset)
            else:
                main_log.info(f'Replacing {short_name} dataset')
                self.current_data[short_name] = dataset

        self.current_data[short_name] = dataset
        # self.save_data(dataset,short_name,folder='main')

        # Process the relaxation data
        if threaded:

            relax_process_worker = Worker(relax_process, dataset)
            relax_process_worker.signals.result.connect(self.post_process_relax)

            self.threadpool.start(relax_process_worker)
        else:
            # pass
            
            self.post_process_relax(relax_process(dataset),test_length=False)

    def post_process_relax(self, fitresult, test_length=True):
        
        seq_name = fitresult.dataset.seq_name
        short_name = short_name_dict[seq_name]
        self.current_results[short_name] = fitresult

        if test_length:
            test_result = fitresult.check_decay()
            test_dt = fitresult.axis[1].values - fitresult.axis[0].values
            test_dt *= 1e3
        else:
            test_result = 0

        if test_result != 0:
            if test_result == -1:  # The trace needs to be longer
                new_dt = epr.round_step(test_dt*2, self.waveform_precision)
                new_tmin = fitresult.axis[-1].values
                new_tmin += new_dt*1e-3
                main_log.info(f"Relaxation trace needs to be longer, setting new dt {new_dt:.2f} ns")
            elif test_result == 1:  # The trace needs to be shorter
                new_dt = epr.round_step(test_dt/2, self.waveform_precision)
                new_tmin = fitresult.axis[0].values
                main_log.info(f"Relaxation trace needs to be shorter, setting new dt {new_dt:.2f} ns")

            nAvgs = fitresult.dataset.attrs['nAvgs']

            if self.worker is not None:
                self.worker.rerun_relax(short_name,dt=new_dt, tmin=new_tmin, averages=nAvgs, autoStop=False, autoIFGain=False)

        else:
            self.save_data(fitresult.dataset,folder='main')
            # self.initialise_deer_settings()
            self.update_deer_settings(update_pulses=False, threaded=False)
            

        if self.waitCondition is not None:  # Wake up the runner thread
                self.waitCondition.wakeAll()
        self.refresh_relax_figure()



    def advanced_mode_inputs(self):
        self.Exp_types.addItems(['auto', '5pDEER', '4pDEER', 'Ref2D'])
        self.bg_model_combo.addItems(list(BackgroundModels.keys()))
        self.ExcPulseSelect.addItems(['Auto', 'Rectangular', 'Gauss'])
        self.RefPulseSelect.addItems(['Auto', 'Rectangular', 'Gauss'])
        self.PumpPulseSelect.addItems(['Auto', 'Rectangular', 'Chirp', 'HS', 'Gauss'])

    def update_quickdeer(self, dataset=None):
        if dataset is None:
            dataset = self.current_data['quickdeer']
        else:
            self.current_data['quickdeer'] = dataset
            self.save_data(dataset,'DEER_init',folder='main')

        self.Tab_widget.setCurrentIndex(3)
        self.q_DEER.current_data['quickdeer'] = dataset
        self.q_DEER.update_inputs_from_dataset()
        # self.q_DEER.update_figure()
        bg_model = BackgroundModels[self.bg_model_combo.currentText()]
        def update_func(x): # This function is called after the DEER analysis is complete
            self.current_results['quickdeer'] = x
            
            if x.MNR < 4:
                main_log.critical(f"QuickDEER MNR is far too low {x.MNR:.2f}, stopping DEER analysis")
                self.worker.stop()
                return None
            elif x.MNR < 10:
                main_log.warning(f"QuickDEER MNR is too low {x.MNR:.2f}, repeating initial DEER with a shorter tau")
                self.self.label_eff /= 4
                # self.initialise_deer_settings()
                self.update_deer_settings(update_pulses=False, threaded=False)
                time = np.min([self.aim_time,self.remaining_time])
                self.worker.repeat_quickdeer()
            else:
                time = None

            self.update_deer_settings(remaining_time=time)
            # self.update_pulses(self._optimise_pulses_in_background()) #reoptimse pulses in background

        self.q_DEER.process_deeranalysis(background_model=bg_model, wait_condition = self.waitCondition, update_func=update_func)

    def update_longdeer(self, dataset=None):
        if dataset is None:
            dataset = self.current_data['longdeer']
        else:
            self.current_data['longdeer'] = dataset
            self.save_data(dataset,'DEER_final',folder='main')

        self.Tab_widget.setCurrentIndex(5)
        self.longDEER.current_data['quickdeer'] = dataset
        self.longDEER.update_inputs_from_dataset()
        self.longDEER.update_figure()
        
        def update_func(x):
            self.current_results['longdeer'] = x
        
        self.longDEER.process_deeranalysis(wait_condition = self.waitCondition,update_func=update_func)

    def update_reptime(self, dataset=None,threaded = True):
        if dataset is None:
            dataset = self.current_data['reptime']
        else:
            self.current_data['reptime'] = dataset
            self.save_data(dataset,'reptime',folder='main')

        # reptime_analysis = epr.ReptimeAnalysis(dataset,dataset.sequence)
        reptime_analysis = epr.ReptimeAnalysis(dataset)
        reptime_analysis.fit()

        if 'ReptimeRecovery' in self.config['autoDEER']:
            ReptimeRecovery = self.config['autoDEER']['ReptimeRecovery']
        else:
            ReptimeRecovery = 0.8
        opt_reptime = reptime_analysis.calc_optimal_reptime(ReptimeRecovery)

        if (opt_reptime*1e-3 > 8) or (opt_reptime*1e-3 < 0.5):
            main_log.warning(f"Reptime optimisation failed. Setting to default for spin system")
            opt_reptime = 3e3

        self.current_results['reptime'] = reptime_analysis
       
        if self.worker is not None:
            self.worker.update_reptime(opt_reptime)
       
        self.SRTBox.setValue(opt_reptime*1e-3)
        main_log.info(f"Reptime {opt_reptime*1e-3:.2g} ms")
        self.update_reptime_figure()
        if self.waitCondition is not None: # Wake up the runner thread
            self.waitCondition.wakeAll()

    def update_reptime_figure(self):

        fig = self.relax_canvas.figure
        self.relax_ax[-2].cla()

        if  not 'reptime' in self.current_results:
            raise ValueError("No reptime analysis found")
        self.current_results['reptime'].plot(axs=self.relax_ax[-2],fig=fig)

    def timeout(self):
        """
        Creates a pop up box as the experiment has timed out
        """
        msg = f'AutoDEER has timedout. The maximum specified time was {self.MaxTime.value()}'
        QMessageBox.about(self,'Warning!', msg)
        main_log.warning(msg)

    def clear_all(self):
        self.current_results = {}
        self.current_data = {}
        self.pulses = {}
        self.userinput = {}
        self.deer_settings = {'ESEEM':None, 'ExpType':'5pDEER'}
        self.q_DEER.current_data = {}
        self.longDEER.current_data = {}

        # Clear the figures
        self.setup_plot_ax.cla()
        self.setup_figure_canvas.draw()
        if isinstance(self.relax_ax, (list,tuple,np.ndarray)):
            for ax in self.relax_ax:
                ax.cla()
        else:
            self.relax_ax.cla()
        self.relax_canvas.draw()
        self.q_DEER.clear_all()
        self.longDEER.clear_all()

    def RunAutoDEER(self, advanced=False):

        self.clear_all()
        self.operating_mode = 'auto'
        self.fixed_tau = None

        self.starttime = time.time()

        if advanced:
            errors = self.validate_advanced_pulse_options()
            if errors:
                QMessageBox.about(self,'ERORR!', 'There are errors in the advanced pulse options:\n' + '\n'.join(errors))
                main_log.error('Could not run autoDEER due to errors in the advanced pulse options:\n' + '\n'.join(errors))
                return None

            
        

        if self.spectromterInterface is None or self.connected is False:
            QMessageBox.about(self,'ERORR!', 'A interface needs to be connected first!')
            main_log.error('Could not run autoDEER. A interface needs to be connected first!')
            return None
        
        # Block the autoDEER buttons
        self.FullyAutoButton.setEnabled(False)
        self.AdvancedAutoButton.setEnabled(False)
        self.resonatorComboBox.setEnabled(False)
        self.fcDoubleSpinBox.setEnabled(False)


        userinput = {}
        userinput['MaxTime'] = self.MaxTime.value()
        userinput['project'] = self.ProjectName.text()
        userinput['sample'] = self.SampleName.text()
        userinput['comment'] = self.commentLineEdit.text()
        userinput['priority'] = self.priotityComboBox.currentText()
        userinput['label_eff'] = self.LabellingEffSpinBox.value()
        userinput['Temp'] = self.TempValue.value()
        userinput['DEER_update_func'] = self.q_DEER.refresh_deer
        userinput['tp'] = self.Min_tp

        self.userinput = userinput

        self.label_eff = self.userinput['label_eff'] / 100


        if self.priotityComboBox.currentText().lower() == 'single':
            self.operating_mode = 'single'
        else:
            self.operating_mode = self.Exp_types.currentText()

        if self.operating_mode == '4pDEER' and not self.config['autoDEER'].get('2D_Dec',True):
            self.operating_mode = '5pDEER'
        
        if advanced:
            self.read_advanced_pulse_options()
            self.read_advanced_sequence_options()
            # Determine if sequence parameters are fixed:
            if (self.AdvSeqOptions['ExpType'].lower() != 'auto') and (self.AdvSequenceCheck.isChecked()):
                
                self.deer_settings['ExpType'] = self.AdvSeqOptions['ExpType']
                if self.AdvSeqOptions['tau1'] > 0 or self.AdvSeqOptions['tau2'] > 0:
                    # Using Fixed taus
                    self.fixed_tau = {}
                    self.fixed_tau['tau1'] = self.AdvSeqOptions['tau1']
                    self.fixed_tau['tau2'] = self.AdvSeqOptions['tau2']
                    
                    if 'tau3' in self.AdvSeqOptions:
                        self.fixed_tau['tau3'] = self.AdvSeqOptions['tau3']
                    self.deer_settings.update(self.fixed_tau)
                    main_log.info(f"tau1 set to {self.deer_settings['tau1']:.2f} us")
                    main_log.info(f"tau2 set to {self.deer_settings['tau2']:.2f} us")
                    self.deer_settings['dt'] = self.AdvSeqOptions['dt']
                else:
                    # Only Setting the sequence
                    self.fixed_tau = None
                
                self.deer_settings['criteria'] = self.priorties[self.userinput['priority']]
                self.deer_settings['autoStop'] = self.Time_autoStop_checkbox.isChecked()
                main_log.info(f"DEER Sequence set to {self.deer_settings['ExpType']}")
                main_log.info(f"tau1 set to {self.deer_settings['tau1']:.2f} us")
                main_log.info(f"tau2 set to {self.deer_settings['tau2']:.2f} us")
                main_log.info(f"DEER Sequence set to {self.deer_settings['ExpType']}")
            
            # else:
            #     self.fixed_tau = None

        
        
        if not 'ESEEM' in self.deer_settings:
            self.deer_settings['ESEEM'] = None

        if self.operating_mode in ['4pDEER','Ref2D']:
            self.create_relax_figure(4)
        else:
            self.create_relax_figure(3)


        try:
            night_hours = self.config['autoDEER']['Night Hours']
        except KeyError:
            night_hours = None

        self.waitCondition = QtCore.QWaitCondition()
        mutex = QtCore.QMutex()
        self.worker = autoDEERWorker(
            self.spectromterInterface, wait=self.waitCondition, mutex=mutex,
            pulses=self.pulses, results=self.current_results, AWG=self.AWG, freq=self.LO, gyro = self.gyro,
            user_inputs=userinput, operating_mode = self.operating_mode, fixed_tau=self.fixed_tau,
            cores=self.cores, night_hours=night_hours,
            )
        

        self.worker.setQ(self.Q)
        

        # self.worker.update_deersettings(self.deer_settings)
        self.worker.signals.update_deer_settings.emit(self.deer_settings)
    
        self.worker.signals.status.connect(self.msgbar.setText)
        self.worker.signals.status.connect(main_log.info)
        self.worker.signals.fsweep_result.connect(self.update_fieldsweep)
        self.worker.signals.fsweep_result.connect(lambda x: self.save_data(x))
        self.worker.signals.respro_result.connect(self.update_respro)
        self.worker.signals.respro_result.connect(lambda x: self.save_data(x))
        # self.worker.signals.optimise_pulses.connect(self.optimise_pulses)
        self.worker.signals.relax_result.connect(self.update_relax)
        self.worker.signals.relax_result.connect(lambda x: self.save_data(x))
        self.worker.signals.T2_result.connect(self.update_relax)
        self.worker.signals.T2_result.connect(lambda x: self.save_data(x))

        self.worker.signals.Relax2D_result.connect(self.update_relax2D)
        self.worker.signals.Relax2D_result.connect(lambda x: self.save_data(x))

        self.worker.signals.quickdeer_result.connect(self.update_quickdeer)
        self.worker.signals.quickdeer_result.connect(lambda x: self.save_data(x))
        self.worker.signals.quickdeer_update.connect(self.q_DEER.refresh_deer)
        self.worker.signals.longdeer_update.connect(self.longDEER.refresh_deer)
        self.worker.signals.longdeer_result.connect(lambda x: self.save_data(x))

        self.worker.signals.longdeer_result.connect(self.update_longdeer)
        self.worker.signals.reptime_scan_result.connect(self.update_reptime)

        self.worker.signals.timeout.connect(self.timeout)
        self.worker.signals.finished.connect(lambda: self.FullyAutoButton.setEnabled(True))
        self.worker.signals.finished.connect(lambda: self.AdvancedAutoButton.setEnabled(True))
        self.worker.signals.finished.connect(lambda: self.resonatorComboBox.setEnabled(True))
        self.worker.signals.finished.connect(lambda: self.fcDoubleSpinBox.setEnabled(True))

        self.stopButton.clicked.connect(self.stopExperiment)

        self.detect_recent_experiments(self.worker)
        time.sleep(2)


        self.threadpool.start(self.worker)
        main_log.info(f"Starting autoDEER")

        return self.worker
    
    def RunFullyAutoDEER(self):
        return self.RunAutoDEER(advanced=False)

    def RunAdvancedAutoDEER(self):
        return self.RunAutoDEER(advanced=True)
    
    def stopExperiment(self):
        
        if hasattr(self, 'spectromterInterface') and self.spectromterInterface is not None:
            self.spectromterInterface.terminate()

        if hasattr(self, 'worker') and self.worker is not None:
            self.worker.stop()

    def create_report(self):
        save_path = QFileDialog.getSaveFileName(self, 'Save File', self.current_folder, ".pdf")
        save_path = save_path[0] + save_path[1]
        report = ad.Reporter(filepath=save_path,pagesize='A4')

        report.add_title('title','autoDEER Report')
        if self.connected:
            report.add_new_section('spec',' Spectrometer')
            report.add_text('spec', 'Local Name: ' + self.config['Spectrometer']['Local Name'])
            report.add_text('spec', 'Manufacturer: ' + self.config['Spectrometer']['Manufacturer'])
            report.add_text('spec', 'Model: ' + self.config['Spectrometer']['Model'])

            report.add_text('spec', 'Resonator: ' + self.resonatorComboBox.currentText())
            report.add_space('spec', height=10)
        
        report.add_new_section('inputs',' User Inputs')
        report.add_text('inputs', 'Sample Name: ' + self.SampleName.text())
        report.add_text('inputs', 'Temperature: ' + str(self.TempValue.value()) + ' K')
        report.add_text('inputs', 'Max Time: ' + str(self.MaxTime.value()) + ' hours')

        report.add_page_break('inputs')

        if 'fieldsweep' in self.current_results:
            report.add_new_section('fieldsweep',' Field Sweep')
            report.add_figure('fieldsweep', self.fsweep_canvas.figure)
            report.add_text('fieldsweep', f"Gyromagnetic Ratio: {self.current_results['fieldsweep'].gyro:.3g} GHz/G")
            if hasattr(self.current_results['fieldsweep'], 'results'):
                report.add_space('fieldsweep')
                report.add_code_block('fieldsweep', self.current_results['fieldsweep'].results.__str__(), title='Fit Results')
            report.add_page_break('fieldsweep')

        if 'respro' in self.current_results:
            report.add_new_section('respro',' Resonator Profile')
            fig,axs = plt.subplots(1,1,figsize=(5, 5))
            fitresult = self.current_results['respro']
            if 'fieldsweep'in self.current_results:
                fitresult.plot(fieldsweep=self.current_results['fieldsweep'],axs=axs, fig=fig);
            else:
                fitresult.plot(axs=axs,fig=fig)


            report.add_figure('respro', fig)
            if hasattr(self.current_results['respro'], 'results'):
                report.add_space('respro')
                report.add_code_block('respro', self.current_results['respro'].results.__str__(), title='Fit Results')
            report.add_page_break('respro')

        if self.pulses != {}:
            pump_pulse = self.pulses['pump_pulse']
            exc_pulse = self.pulses['exc_pulse']
            ref_pulse = self.pulses['ref_pulse']
            report.add_new_section('pulses',' Optimised Pulses')
            fig,axs = plt.subplots(1,1,figsize=(5, 5))
            ad.plot_overlap(self.current_results['fieldsweep'], pump_pulse, exc_pulse,ref_pulse, axs=axs,fig=fig)

            report.add_figure('pulses', fig)
            report.add_space('pulses')

            report.add_code_block('pulses', exc_pulse.__str__(), title='Excitation Pulse')
            report.add_code_block('pulses', ref_pulse.__str__(), title='Refocusing Pulse')
            report.add_code_block('pulses', pump_pulse.__str__(), title='Pump Pulse')


        if 'CP-relax' in self.current_results:
            report.add_new_section('relax',' Relaxation')
            report.add_figure('relax', self.relax_canvas.figure)
            if hasattr(self.current_results['CP-relax'], 'results'):
                report.add_space('relax')
                report.add_code_block('relax', self.current_results['CP-relax'].results.__str__(), title='Fit Results')
            report.add_page_break('relax') 


        if 'quickdeer' in self.q_DEER.current_data:
            report.add_new_section('quickdeer',' QuickDEER')
            fig,axs = plt.subplot_mosaic([['Primary_time'], 
                                          ['Primary_dist']], figsize=(6,6))
            
            ad.DEERanalysis_plot(self.q_DEER.fitresult, background=True, ROI=self.q_DEER.fitresult.ROI, axs= axs, fig=fig,text=False)
            report.add_figure('quickdeer', fig)
            
            if 'quickdeer' in self.current_results:
                report.add_space('quickdeer')
                report.add_code_block('quickdeer', self.current_results['quickdeer'].__str__(), title='Fit Results')
            report.add_page_break('quickdeer')
            
        report._build()
        pass

def main():


    app = QApplication([])
    app.setWindowIcon(QtGui.QIcon('icons:Square_logo.png'))
    app.setApplicationName('autoDEER')
    app.setApplicationDisplayName('autoDEER')
    app.setApplicationVersion(str(ad.__version__))
    apply_stylesheet(app, theme='light_purple.xml')
    app.setDesktopFileName('autoDEER')
    window = autoDEERUI()
    window.show()
    app.exec()

if __name__ == '__main__':
    main()
    
