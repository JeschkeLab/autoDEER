from PyQt6.QtWidgets import QApplication, QMainWindow, QFileDialog,QMessageBox
from PyQt6 import uic
import PyQt6.QtCore as QtCore
import PyQt6.QtGui as QtGui
from pathlib import Path
import sys, traceback, os
from threadpoolctl import threadpool_limits


from matplotlib.backends.backend_qtagg import FigureCanvas, NavigationToolbar2QT
import matplotlib.pyplot as plt
import autodeer as ad
import numpy as np
from autodeer.gui.tools import *
from autodeer.gui.autoDEER_worker import autoDEERWorker
from autodeer.gui.quickdeer import DEERplot
from autodeer.gui.log import LogDialog
import yaml
import time
import datetime
import logging
from autodeer.Logging import setup_logs, change_log_level

# main_log = logging.getLogger('autoDEER')
from queue import Queue

package_directory = os.path.dirname(os.path.abspath(__file__))

QtCore.QDir.addSearchPath('icons', f"{package_directory}/resources")

SampleConcComboBox_opts = {'Normal': 1, 'High (0.5x)':0.5, 'Low (2x)':2, 'Very Low (5x)':5}

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

def fieldsweep_fit(fsweep_analysis):
    fsweep_analysis.fit(xtol=1e-5, lin_maxiter=100)

    return fsweep_analysis

def respro_process(dataset, fieldsweep=None,cores=1):
    respro = ad.ResonatorProfileAnalysis(
        dataset
    )    
    with threadpool_limits(limits=cores, user_api='blas'):
        respro.fit(cores=cores)

    if fieldsweep is not None:
        LO_new = fieldsweep.LO + ad.optimise_spectra_position(respro, fieldsweep)
        return respro, LO_new


    return respro

def relax_process(dataset):

    # if dataset.axes[0].max()>500:
    #     dataset.axes = [dataset.axes[0]/1e3]
    # else:
    #     dataset.axes = [dataset.axes[0]]
    if dataset['tau1'].max() > 1e4:
        dataset['tau1'] /= 1e3
    CP_data = ad.CarrPurcellAnalysis(dataset)
    CP_data.fit('double')

    return CP_data

def T2_process(dataset):
    CP_data = ad.CarrPurcellAnalysis(dataset)
    CP_data.fit('double')

    return CP_data

class autoDEERUI(QMainWindow):


    def __init__(self):
        super().__init__()
 
        # loading the ui file with uic module
        uic.loadUi(f"{package_directory}/gui2.ui", self)
        logo_pixmap = QtGui.QPixmap('icons:logo.png')
        logo_pixmap = logo_pixmap.scaledToHeight(60)
        self.logo.setPixmap(logo_pixmap)
        self.set_spectrometer_connected_light(0)

        self.setWindowTitle("autoDEER")


        self.qDEER_tab.layout().addWidget(DEERplot())
        self.q_DEER = self.qDEER_tab.layout().itemAt(0).widget()
        self.lDEER_tab.layout().addWidget(DEERplot())
        self.longDEER = self.lDEER_tab.layout().itemAt(0).widget()


        self.fsweep_toolbar()
        self.respro_toolbar()
        self.relax_toolbar()
        self.create_respro_figure()
        self.create_relax_figure()
        self.create_fieldsweep_figure()
        self.advanced_mode_inputs()

        self.threadpool = QtCore.QThreadPool()
        self.current_results = {}
        self.current_data = {}
        self.pulses = {}
        self.spectromterInterface = None
        self.waitCondition = None
        self.queue = Queue()

        self.FullyAutoButton.clicked.connect(self.RunFullyAutoDEER)
        self.AdvancedAutoButton.clicked.connect(self.RunAdvanedAutoDEER)

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
        self.show_respro.clicked.connect(lambda: self.update_respro())
        self.OptimisePulsesButton.clicked.connect(lambda: self.optimise_pulses_button())

        self.SampleConcComboBox.addItems(list(SampleConcComboBox_opts.keys()))
        

        self.Tab_widget.setCurrentIndex(0)
        # set current folder to home directory
        self.current_folder = None
        self.config = None
        self.connected = False
        self.DL_params = {}
        

        self.LO = 0
        self.gyro = 0.002803632236095
        self.cores = 1

    def set_spectrometer_connected_light(self, state):
        if state == 0:
            light_pixmap = QtGui.QPixmap('icons:Red.png')
        elif state == 1:
            light_pixmap = QtGui.QPixmap('icons:Green.png')
        elif state == 2:
            light_pixmap = QtGui.QPixmap('icons:Yellow.png')
        
        light_pixmap = light_pixmap.scaledToHeight(30)
        self.Connected_Light.setPixmap(light_pixmap)

    def load_folder(self):
        file = str(QFileDialog.getExistingDirectory(self, "Select Directory",str(Path.home())))
        self.pathLineEdit.setText(file)
        self.current_folder = file

        setup_logs(self.current_folder)
        global main_log
        main_log = logging.getLogger('autoDEER')
        main_log.info(f"Loading folder {file}")
    
    def load_epr_file(self, store_location):

        filename, _= QFileDialog.getOpenFileName(
            self,"Select a File", self.current_folder,"Data (*.DTA *.mat)")
        
        if filename:
                path = Path(filename)
                filename_edit = str(path)

        dataset = ad.eprload(filename_edit)
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
                from autodeer.hardware.dummy import dummyInterface
                self.spectromterInterface = dummyInterface(filename_edit)
            elif model == 'ETH_AWG':
                from autodeer.hardware.ETH_awg import ETH_awg_interface
                self.spectromterInterface = ETH_awg_interface()
            elif model == 'Bruker_MPFU':
                from autodeer.hardware.Bruker_MPFU import BrukerMPFU
                self.spectromterInterface = BrukerMPFU(filename_edit)
                self.spectromterInterface.savefolder = self.current_folder
            elif model == 'Bruker_AWG':
                from autodeer.hardware.Bruker_AWG import BrukerAWG
                self.spectromterInterface = BrukerAWG(filename_edit)
        except ImportError:
            QMessageBox.about(self,'ERORR!', 
                              'The spectrometer interface could not be loaded!\n'+
                              'Please check that the correct packages are installed!\n'+
                              'See the documentation for more information.')
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

        # Get user preferences
        try:
            self.cores = int(self.config['autoDEER']['cores'])
            self.q_DEER.cores = self.cores
            self.longDEER.cores = self.cores
        except:
            self.q_DEER.cores = 1
            self.longDEER.cores = 1

        # Extract DeerLab fit parameters
        if ('autoDEER' in self.config) and ('DeerLab' in self.config['autoDEER']):
            self.DL_params = self.config['autoDEER']['DeerLab']
        else:
            self.DL_params = {}

        self.q_DEER.DL_params = self.DL_params
        self.longDEER.DL_params = self.DL_params

        try:
            self.waveform_precision = spectrometer['Bridge']['Pulse dt']
        except KeyError:
            self.waveform_precision = 1 


    def select_resonator(self):
        key = self.resonatorComboBox.currentText()
        main_log.info(f"Selecting resonator {key}")
        self.LO = self.config['Resonators'][key]['Center Freq']
        self.fcDoubleSpinBox.setValue(self.LO)
        main_log.info(f"Setting LO to {self.LO} GHz")
    
    def change_LO(self):
        self.LO = self.fcDoubleSpinBox.value()
        main_log.info(f"Setting LO to {self.LO} GHz")

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

    def save_data(self,dataset,filename):
        timestamp = datetime.datetime.now().strftime(r'%Y%m%d_%H%M_')
        fullname = timestamp + self.SampleName.text() + '_' +filename +'.h5'
        dataset.to_netcdf(os.path.join(self.current_folder,fullname),engine='h5netcdf',invalid_netcdf=True)
        main_log.debug(f"Saving dataset to {fullname}")

    def fsweep_toolbar(self):
        upload_icon = QtGui.QIcon('icons:upload.png')
        self.Load_button.setIcon(upload_icon)
        self.Load_button.clicked.connect(lambda x: self.load_epr_file('fieldsweep'))

        refresh_icon = QtGui.QIcon('icons:refresh.png')
        self.Refresh_button.setIcon(refresh_icon)
        self.Refresh_button.clicked.connect(lambda: self.update_fieldsweep())

    def respro_toolbar(self):
        upload_icon = QtGui.QIcon('icons:upload.png')
        self.respro_Load_button.setIcon(upload_icon)
        self.respro_Load_button.clicked.connect(lambda: self.load_epr_file('respro'))

        refresh_icon = QtGui.QIcon('icons:refresh.png')
        self.respro_Refresh_button.setIcon(refresh_icon)
        self.respro_Refresh_button.clicked.connect(lambda: self.update_respro())

    def relax_toolbar(self):
        upload_icon = QtGui.QIcon('icons:upload.png')
        self.relax_Load_button.setIcon(upload_icon)
        self.relax_Load_button.clicked.connect(lambda: self.load_epr_file('relax'))

        refresh_icon = QtGui.QIcon('icons:refresh.png')
        self.relax_Refresh_button.setIcon(refresh_icon)
        self.relax_Refresh_button.clicked.connect(lambda: self.update_relax())

    def update_fieldsweep(self, dataset=None):

        if dataset is None:
            dataset = self.current_data['fieldsweep']
        else:
            self.current_data['fieldsweep'] = dataset

        fsweep_analysis = ad.FieldSweepAnalysis(dataset)
        fsweep_analysis.calc_gyro()
        main_log.info(f"Calculated gyro {fsweep_analysis.gyro*1e3:.3f} MHz/T")
        self.fsweep_ax.cla()
        fsweep_analysis.plot(axs=self.fsweep_ax,fig=self.fsweep_canvas.figure)
        self.fsweep_canvas.draw()
        self.Tab_widget.setCurrentIndex(1)

        worker = Worker(fieldsweep_fit, fsweep_analysis)
        worker.signals.result.connect(self.refresh_fieldsweep_after_fit)
        
        self.threadpool.start(worker)

    def create_fieldsweep_figure(self):
        fig, axs  = plt.subplots(1,1,figsize=(12.5, 6.28))
        self.fsweep_canvas = FigureCanvas(fig)
        self.fsweep_v_left.addWidget(self.fsweep_canvas)
        Navbar = NavigationToolbar2QT(self.fsweep_canvas, self)
        Navbar.setMaximumHeight(24)
        self.fsweep_v_left.addWidget(Navbar)
        self.fsweep_ax = axs
 
    def refresh_fieldsweep_after_fit(self, fitresult):

        self.current_results['fieldsweep'] = fitresult
        self.gyro = fitresult.gyro
        if self.waitCondition is not None: # Wake up the runner thread
            self.waitCondition.wakeAll()
        
        self.fsweep_ax.cla()
        fitresult.plot(axs=self.fsweep_ax,fig=self.fsweep_canvas.figure)
        self.fsweep_canvas.draw()


        # Add fit resultss

        self.gyroSpinBox.setValue(fitresult.gyro*1e3)
        self.gxSpinBox.setValue(-0.0025 * fitresult.results.az + 2.0175)
        self.gySpinBox.setValue(fitresult.results.gy)
        self.gzSpinBox.setValue(fitresult.results.gz)
        self.AxSpinBox.setValue(fitresult.results.axy*28.0328)
        self.AySpinBox.setValue(fitresult.results.axy*28.0328)
        self.AzSpinBox.setValue(fitresult.results.az)
        self.GBSpinBox.setValue(fitresult.results.GB)
        self.BoffsetSpinBox.setValue(fitresult.results.Boffset)

        gxCI = -0.0025 *fitresult.results.azUncert.ci(95)+ 2.0175
        self.gxCI.setText(f"({gxCI[0]:.4f},{gxCI[1]:.4f})")
        self.gyCI.setText(getCIstring(fitresult.results.gyUncert))
        self.gzCI.setText(getCIstring(fitresult.results.gzUncert))
        self.Tab_widget.setCurrentIndex(1)



    def update_respro(self, dataset=None):

        if dataset is None:
            dataset = self.current_data['respro']
        else:
            self.current_data['respro'] = dataset

        if hasattr(dataset, 'sequence'):
            f_axis = dataset.sequence.LO.value + dataset.sequence.LO.axis[0]['axis']# - 0.3 # Fix this necessary offset
        else:
            # Assuming mat file
            # TODO: Change to new data format
            f_axis = dataset['LO'].data #- 0.3
            # f_axis = dataset.params['parvars'][2]['vec'][:,1] + dataset.params['LO']


        # worker = Worker(respro_process, dataset, f_axis,self.current_results['fieldsweep'], cores=self.cores)
        worker = Worker(respro_process, dataset, self.current_results['fieldsweep'], cores=self.cores)

        worker.signals.result.connect(self.refresh_respro)

        self.threadpool.start(worker)

    def create_respro_figure(self):
        fig, axs  = plt.subplots(1,1,figsize=(12.5, 6.28))
        self.respro_canvas = FigureCanvas(fig)
        Navbar = NavigationToolbar2QT(self.respro_canvas, self)
        Navbar.setMaximumHeight(24)
        self.respro_v_left.addWidget(self.respro_canvas)
        self.respro_v_left.addWidget(Navbar)
        self.respro_ax = axs


    def refresh_respro(self, *args):

        

        if len(args[0]) == 2:
            fitresult = args[0][0]
            self.LO = args[0][1]
            if self.worker is not None:
                self.worker.update_LO(self.LO)
            print(f"New LO frequency: {self.LO:.2f} GHz")
            main_log.info(f"Setting LO to {self.LO:.6f} GHz")
        else:
            fitresult = args[0]

        self.current_results['respro'] = fitresult
        if self.waitCondition is not None: # Wake up the runner thread
            self.waitCondition.wakeAll()
            
        self.respro_ax.cla()

        if 'fieldsweep'in self.current_results:
            fitresult.plot(fieldsweep=self.current_results['fieldsweep'],axs=self.respro_ax,fig=self.respro_canvas.figure);
        else:
            fitresult.plot(axs=self.respro_ax,fig=self.respro_canvas.figure)

        self.respro_canvas.draw()
        # Add fit results

        self.centreFrequencyDoubleSpinBox.setValue(fitresult.results.fc)
        self.centreFrequencyCI.setText(f"({fitresult.results.fcUncert.ci(95)[0]:.2f},{fitresult.results.fcUncert.ci(95)[1]:.2f})")
        self.qDoubleSpinBox.setValue(fitresult.results.q)
        self.qCI.setText(f"({fitresult.results.qUncert.ci(95)[0]:.2f},{fitresult.results.qUncert.ci(95)[1]:.2f})")
        self.Tab_widget.setCurrentIndex(2)
        main_log.info(f"Resonator centre frequency {fitresult.results.fc:.4f} GHz")
        self.optimise_pulses_button()


    def optimise_pulses_button(self):
        if (self.pulses is None) or (self.pulses == {}):
            self.optimise_pulses()
        elif self.pulses['pump_pulse'] is None:
            self.optimise_pulses()
        else:
            self.update_optimise_pulses_figure()
        

    def optimise_pulses(self, pulses=None):
        if (pulses is None) or pulses == {}:
            self.pulses = ad.build_default_pulses(self.AWG)
            # pump_pulse = ad.HSPulse(tp=120, init_freq=-0.25, final_freq=-0.03, flipangle=np.pi, scale=0, order1=6, order2=1, beta=10)
            # exc_pulse = ad.RectPulse(tp=16, freq=0.02, flipangle=np.pi/2, scale=0)
            # ref_pulse = exc_pulse.copy(flipangle=np.pi)
        
        pump_pulse = self.pulses['pump_pulse']
        ref_pulse = self.pulses['ref_pulse']
        exc_pulse = self.pulses['exc_pulse']

        pump_pulse, exc_pulse, ref_pulse = ad.optimise_pulses(self.current_results['fieldsweep'], pump_pulse, exc_pulse, ref_pulse)
        
        # self.pulses = {'pump_pulse':pump_pulse, 'exc_pulse':exc_pulse, 'ref_pulse':ref_pulse}    
        self.pulses['pump_pulse'] = pump_pulse
        self.pulses['ref_pulse'] = ref_pulse
        self.pulses['exc_pulse'] = exc_pulse
        self.pulses['det_event'].freq = exc_pulse.freq
        main_log.info(f"Optimised pulses")

        self.worker.new_pulses(self.pulses)
        
        if self.waitCondition is not None: # Wake up the runner thread
            self.waitCondition.wakeAll()

        self.update_optimise_pulses_figure()
        self.Tab_widget.setCurrentIndex(2)
    
    def update_optimise_pulses_figure(self):
        pump_pulse = self.pulses['pump_pulse']
        exc_pulse = self.pulses['exc_pulse']
        ref_pulse = self.pulses['ref_pulse']

        self.respro_ax.cla()
        ad.plot_overlap(self.current_results['fieldsweep'], pump_pulse, exc_pulse,ref_pulse, axs=self.respro_ax,fig=self.respro_canvas.figure)
        self.respro_canvas.draw()

        # update the pulse parameter grid
        type_to_pulse_hash = {ad.RectPulse:'Rect', ad.ChirpPulse:'Chirp', ad.HSPulse:'HS'}
        if isinstance(exc_pulse, ad.RectPulse):
            self.ExcFreqBox.setValue(param_in_MHz(exc_pulse.freq))
            self.ExcBWBox.setValue(exc_pulse.tp.value)
            self.ExcBWBox.setSuffix(' ns')
            self.ExcTypeLine.setText('Rect')
        else:
            center_freq = (param_in_MHz(exc_pulse.final_freq) + param_in_MHz(exc_pulse.init_freq))/2
            self.ExcFreqBox.setValue(center_freq)
            self.ExcBWBox.setValue(param_in_MHz(exc_pulse.BW))
            self.ExcBWBox.setSuffix(' MHz')
            self.ExcTypeLine.setText(type_to_pulse_hash[type(exc_pulse)])
        
        if isinstance(ref_pulse, ad.RectPulse):
            self.RefFreqBox.setValue(param_in_MHz(ref_pulse.freq))
            self.RefBWBox.setValue(ref_pulse.tp.value)
            self.RefBWBox.setSuffix(' ns')
            self.RefTypeLine.setText('Rect')
        else:
            center_freq = (param_in_MHz(ref_pulse.final_freq) + param_in_MHz(ref_pulse.init_freq))/2
            self.RefFreqBox.setValue(center_freq)
            self.RefBWBox.setValue(param_in_MHz(ref_pulse.BW))
            self.RefBWBox.setSuffix(' MHz')
            self.RefTypeLine.setText(type_to_pulse_hash[type(ref_pulse)])
        
        if isinstance(pump_pulse, ad.RectPulse):
            self.PumpFreqBox.setValue(param_in_MHz(pump_pulse.freq))
            self.PumpBWBox.setValue(pump_pulse.tp.value)
            self.PumpBWBox.setSuffix(' ns')
            self.PumpTypeLine.setText('Rect')
        else:
            center_freq = (param_in_MHz(pump_pulse.final_freq) + param_in_MHz(pump_pulse.init_freq))/2
            self.PumpFreqBox.setValue(center_freq)
            self.PumpBWBox.setValue(param_in_MHz(pump_pulse.BW))
            self.PumpBWBox.setSuffix(' MHz')
            self.PumpTypeLine.setText(type_to_pulse_hash[type(pump_pulse)])

    def update_relax(self, dataset=None):
        if dataset is None:
            dataset = self.current_data['relax']
        else:
            self.current_data['relax'] = dataset

        worker = Worker(relax_process, dataset)
        worker.signals.result.connect(self.refresh_relax)

        self.threadpool.start(worker)

    def create_relax_figure(self):
        fig, axs  = plt.subplots(1,1,figsize=(12.5, 6.28))
        self.relax_canvas = FigureCanvas(fig)
        Navbar = NavigationToolbar2QT(self.relax_canvas, self)
        Navbar.setMaximumHeight(24)
        self.relax_v_left.addWidget(self.relax_canvas)
        self.relax_v_left.addWidget(Navbar)
        self.relax_ax = axs

    def refresh_relax_figure(self):
        self.relax_ax.cla()

        fig = self.relax_canvas.figure
        axs = self.relax_ax

        axs.set_xlabel("Time ($\mu s$)")
        axs.set_ylabel("Signal (A. U.)")
        
        if 'relax' in self.current_results:
            CP_results = self.current_results['relax']
            CP_data = CP_results.data
            CP_data /= CP_data.max()
            if hasattr(CP_results, "fit_result"):
                axs.plot(CP_results.axis, CP_data, '.', label='CP', color='C1', ms=6)
                axs.plot(CP_results.axis, CP_results.func(
                    CP_results.axis, *CP_results.fit_result[0]), label='CP-fit', color='C1', lw=2)
            else:
                axs.plot(CP_results.axis, CP_data, label='data', color='C1')
        
        if 'T2_relax' in self.current_results:
            T2_results = self.current_results['T2_relax']
            T2_data = T2_results.data
            T2_data /= T2_data.max()
            if hasattr(T2_results, "fit_result"):
                axs.plot(T2_results.axis, T2_data, '.', label='T2', color='C2', ms=6)
                axs.plot(T2_results.axis, T2_results.func(
                    T2_results.axis, *T2_results.fit_result[0]), label='T2-fit', color='C2', lw=2)
            else:
                axs.plot(T2_results.axis, T2_data, label='data', color='C2')
        self.relax_canvas.draw()

    def refresh_relax(self, fitresult):
        self.current_results['relax'] = fitresult

        # self.relax_ax.cla()
        # fitresult.plot(axs=self.relax_ax,fig=self.relax_canvas.figure)
        # self.relax_canvas.draw()

        self.refresh_relax_figure()

        tau2hrs = fitresult.find_optimal(SNR_target=45/0.5, target_time=2, target_step=0.015)
        max_tau = fitresult.find_optimal(SNR_target=45/0.5, target_time=24, target_step=0.015)
        self.current_results['relax'].tau2hrs = tau2hrs
        self.current_results['relax'].max_tau = max_tau
        self.DipolarEvoMax.setValue(max_tau)
        self.DipolarEvo2hrs.setValue(tau2hrs)
        self.Tab_widget.setCurrentIndex(3)
        main_log.info(f"2hr dipolar evo {tau2hrs:.2f} us")
        main_log.info(f"Max dipolar evo {max_tau:.2f} us")

        if self.waitCondition is not None: # Wake up the runner thread
            self.waitCondition.wakeAll()
        
    
    def update_T2(self,dataset=None):
        if dataset is None:
            dataset = self.current_data['T2_relax']
        else:
            self.current_data['T2_relax'] = dataset
        
        # Since the T2 values are not used for anything there is no need pausing the spectrometer    
        if self.waitCondition is not None: # Wake up the runner thread
            self.waitCondition.wakeAll()

        worker = Worker(T2_process, dataset)
        worker.signals.result.connect(self.refresh_T2)

        self.threadpool.start(worker)

    def refresh_T2(self, fitresult):
        self.current_results['T2_relax'] = fitresult
        self.refresh_relax_figure()

        
    def advanced_mode_inputs(self):
        self.Exp_types.addItems(['5pDEER','4pDEER','3pDEER','nDEER'])
        self.ExcPulseSelect.addItems(['Auto', 'Rectangular','Chirp','HS', 'Gauss'])
        self.RefPulseSelect.addItems(['Auto', 'Rectangular','Chirp','HS', 'Gauss'])
        self.PumpPulseSelect.addItems(['Auto', 'Rectangular','Chirp','HS', 'Gauss'])

    def update_quickdeer(self, dataset=None):
        if dataset is None:
            dataset = self.current_data['quickdeer']
        else:
            self.current_data['quickdeer'] = dataset

        self.Tab_widget.setCurrentIndex(4)
        self.q_DEER.current_data['quickdeer'] = dataset
        self.q_DEER.update_inputs_from_dataset()
        self.q_DEER.update_figure()
        def update_func(x):
            self.current_results['quickdeer'] = x
            rec_tau = self.current_data['quickdeer'].rec_tau_max
            dt = self.current_data['quickdeer'].rec_dt * 1e3
            dt = ad.round_step(dt,self.waveform_precision)
            max_tau = self.results['relax'].max_tau
            tau = np.min([rec_tau,max_tau])
            tau1 = ad.round_step(tau,self.waveform_precision)
            tau2 = ad.round_step(tau,self.waveform_precision)
            tau3 = 0.3
            self.worker.update_deersettings(
                tau1=tau1,
                tau2=tau2,
                tau3=tau3,
                dt=dt
            )
        self.q_DEER.process_deeranalysis(wait_condition = self.waitCondition,update_func=update_func)

    def update_longdeer(self, dataset=None):
        if dataset is None:
            dataset = self.current_data['longdeer']
        else:
            self.current_data['longdeer'] = dataset

        self.Tab_widget.setCurrentIndex(4)
        self.longDEER.current_data['quickdeer'] = dataset
        self.longDEER.update_inputs_from_dataset()
        self.longDEER.update_figure()
        
        def update_func(x):
            self.current_results['longdeer'] = x
        
        self.longDEER.process_deeranalysis(wait_condition = self.waitCondition,update_func=update_func)

    def update_reptime(self, dataset=None):
        if dataset is None:
            dataset = self.current_data['reptime']
        else:
            self.current_data['reptime'] = dataset

        # reptime_analysis = ad.ReptimeAnalysis(dataset,dataset.sequence)
        reptime_analysis = ad.ReptimeAnalysis(dataset)
        reptime_analysis.fit()
        opt_reptime = reptime_analysis.calc_optimal_reptime(0.8)

        self.current_results['reptime'] = reptime_analysis
        if self.worker is not None:
            self.worker.update_reptime(opt_reptime)
        main_log.info(f"Reptime {opt_reptime*1e-3:.2g} ms")
        if self.waitCondition is not None: # Wake up the runner thread
            self.waitCondition.wakeAll()
    
    def timeout(self):
        """
        Creates a pop up box as the experiment has timed out
        """
        msg = f'AutoDEER has timedout. The maximum specified time was {self.MaxTime.value()}'
        QMessageBox.about(self,'Warning!', msg)
        main_log.warning(msg)

    def RunFullyAutoDEER(self):

        if self.spectromterInterface is None or self.connected is False:
            QMessageBox.about(self,'ERORR!', 'A interface needs to be connected first!')
            main_log.error('Could not run autoDEER. A interface needs to be connected first!')
            return None

        userinput = {}
        userinput['MaxTime'] = self.MaxTime.value()
        userinput['sample'] = self.SampleName.text()
        userinput['Temp'] = self.TempValue.value()
        userinput['DEER_update_func'] = self.q_DEER.refresh_deer
        userinput['SampleConc'] = SampleConcComboBox_opts[self.SampleConcComboBox.currentText()]

        # Block the autoDEER buttons
        self.FullyAutoButton.setEnabled(False)
        self.AdvancedAutoButton.setEnabled(False)
        self.resonatorComboBox.setEnabled(False)
        self.fcDoubleSpinBox.setEnabled(False)

        self.waitCondition = QtCore.QWaitCondition()
        mutex = QtCore.QMutex()
        self.worker = autoDEERWorker(
            self.spectromterInterface,wait=self.waitCondition,mutex=mutex,
            pulses=self.pulses,results=self.current_results, AWG=self.AWG, LO=self.LO, gyro = self.gyro,
            user_inputs=userinput, cores=self.cores )
    
        self.worker.signals.status.connect(self.msgbar.setText)
        self.worker.signals.status.connect(main_log.info)
        self.worker.signals.fsweep_result.connect(self.update_fieldsweep)
        self.worker.signals.fsweep_result.connect(lambda x: self.save_data(x,'EDFS'))
        self.worker.signals.respro_result.connect(self.update_respro)
        self.worker.signals.respro_result.connect(lambda x: self.save_data(x,'ResPro'))
        # self.worker.signals.optimise_pulses.connect(self.optimise_pulses)
        self.worker.signals.relax_result.connect(self.update_relax)
        self.worker.signals.relax_result.connect(lambda x: self.save_data(x,'CP_Q'))
        self.worker.signals.T2_result.connect(self.update_T2)
        self.worker.signals.T2_result.connect(lambda x: self.save_data(x,'T2'))
        self.worker.signals.quickdeer_result.connect(self.update_quickdeer)
        self.worker.signals.quickdeer_result.connect(lambda x: self.save_data(x,'DEER_5P_Q_quick'))
        self.worker.signals.quickdeer_update.connect(self.q_DEER.refresh_deer)
        self.worker.signals.longdeer_update.connect(self.longDEER.refresh_deer)
        self.worker.signals.longdeer_result.connect(lambda x: self.save_data(x,'DEER_5P_Q_long'))

        self.worker.signals.longdeer_result.connect(self.update_longdeer)


        self.worker.signals.reptime_scan_result.connect(self.update_reptime)

        self.worker.signals.timeout.connect(self.timeout)
        self.worker.signals.finished.connect(lambda: self.FullyAutoButton.setEnabled(True))
        self.worker.signals.finished.connect(lambda: self.AdvancedAutoButton.setEnabled(True))
        self.worker.signals.finished.connect(lambda: self.resonatorComboBox.setEnabled(True))
        self.worker.signals.finished.connect(lambda: self.fcDoubleSpinBox.setEnabled(True))

        self.stopButton.clicked.connect(self.worker.stop)


        self.threadpool.start(self.worker)
        main_log.info(f"Starting autoDEER")

        return self.worker

    def RunAdvanedAutoDEER(self):

        if self.spectromterInterface is None or self.connected is False:
            QMessageBox.about(self,'ERORR!', 'A interface needs to be connected first!')
            main_log.error('Could not run autoDEER. A interface needs to be connected first!')
            return None
        
        if self.current_folder is None:
            QMessageBox.about(self,'ERORR!', 'A folder needs to be selected first!')
            main_log.error('Could not run autoDEER. A folder needs to be selected first!')
            return None

        
        
        userinput = {}
        userinput['MaxTime'] = self.MaxTime.value()
        userinput['sample'] = self.SampleName.text()
        userinput['Temp'] = self.TempValue.value()
        userinput['ExpType'] = self.Exp_types.currentText()
        userinput['tau1'] = self.Tau1Value.value()
        userinput['tau2'] = self.Tau2Value.value()
        userinput['tau3'] = self.Tau3Value.value()
        userinput['ExcPulse'] = self.ExcPulseSelect.currentText()
        userinput['RefPulse'] = self.RefPulseSelect.currentText()
        userinput['PumpPulse'] = self.PumpPulseSelect.currentText()
        userinput['DEER_update_func'] = self.q_DEER.refresh_deer

        # Block the autoDEER buttons
        self.FullyAutoButton.setEnabled(False)
        self.AdvancedAutoButton.setEnabled(False)

        self.waitCondition = QtCore.QWaitCondition()
        mutex = QtCore.QMutex()


        worker = autoDEERWorker(
            self.spectromterInterface,wait=self.waitCondition,mutex=mutex,
            results=self.current_results,LO=self.LO, gyro = self.gyro,
            user_inputs=userinput, cores=self.cores)
        worker.signals.status.connect(self.msgbar.setText)
        worker.signals.fsweep_result.connect(self.update_fieldsweep)
        worker.signals.respro_result.connect(self.update_respro)
        worker.signals.optimise_pulses.connect(self.optimise_pulses)
        worker.signals.relax_result.connect(self.update_relax)
        worker.signals.quickdeer_result.connect(self.update_quickdeer)
        worker.signals.quickdeer_update.connect(self.q_DEER.refresh_deer)
        worker.signals.longdeer_update.connect(self.longDEER.refresh_deer)

        worker.signals.reptime_scan_result.connect(self.update_reptime)

        worker.signals.finished.connect(lambda: self.FullyAutoButton.setEnabled(True))
        worker.signals.finished.connect(lambda: self.AdvancedAutoButton.setEnabled(True))

        self.worker = worker
        self.threadpool.start(self.worker)
        main_log.info(f"Starting autoDEER with user inputs")

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


        if 'relax' in self.current_results:
            report.add_new_section('relax',' Relaxation')
            report.add_figure('relax', self.relax_canvas.figure)
            if hasattr(self.current_results['relax'], 'results'):
                report.add_space('relax')
                report.add_code_block('relax', self.current_results['relax'].results.__str__(), title='Fit Results')
            report.add_page_break('relax') 


        if 'quickdeer' in self.q_DEER.current_data:
            report.add_new_section('quickdeer',' QuickDEER')
            fig,axs = plt.subplot_mosaic([['Primary_time'], 
                                          ['Primary_dist']], figsize=(6,6))
            
            ad.DEERanalysis_plot(self.q_DEER.fitresult, background=False, ROI=self.q_DEER.fitresult.ROI, axs= axs, fig=fig,text=False)
            report.add_figure('quickdeer', fig)
            
            if 'quickdeer' in self.current_results:
                report.add_space('quickdeer')
                report.add_code_block('quickdeer', self.current_results['quickdeer'].__str__(), title='Fit Results')
            report.add_page_break('quickdeer')
            
        report._build()
        pass

if __name__ == '__main__':
    app = QApplication([])
    app.setWindowIcon(QtGui.QIcon('icons:Square_logo.png'))
    app.setApplicationName('autoDEER')
    app.setApplicationDisplayName('autoDEER')
    app.setApplicationVersion(str(ad.__version__))
    window = autoDEERUI()
    window.show()
    app.exec()
