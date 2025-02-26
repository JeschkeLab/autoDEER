from PyQt6.QtWidgets import QApplication, QMainWindow, QFileDialog,QMessageBox, QDialog, QPushButton,QVBoxLayout
from PyQt6 import uic
import PyQt6.QtCore as QtCore
import PyQt6.QtGui as QtGui
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
import yaml
import time
import datetime
import logging
from autodeer.Logging import setup_logs, change_log_level
from deerlab import store_pickle
import deerlab as dl

main_log = logging.getLogger('autoDEER')
from queue import Queue

package_directory = os.path.dirname(os.path.abspath(__file__))

QtCore.QDir.addSearchPath('icons', f"{package_directory}/resources")

SampleConcComboBox_opts = {'Normal': 1, 'High (0.5x)':0.5, 'Low (2x)':2, 'Very Low (5x)':5}
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

    # if dataset.axes[0].max()>500:
    #     dataset.axes = [dataset.axes[0]/1e3]
    # else:
    #     dataset.axes = [dataset.axes[0]]
    if dataset['tau1'].max() > 1e4:
        dataset['tau1'] /= 1e3
    CP_data = epr.CarrPurcellAnalysis(dataset)
    CP_data.fit('auto')

    return CP_data

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
        self.advanced_mode_inputs()

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

        self.SampleConcComboBox.addItems(list(SampleConcComboBox_opts.keys()))
        

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
        self.gyro = 0.002803632236095
        self.gyro = 0.002808859721083
        self.cores = 1
        self.Min_tp=12

        self.deer_settings = {'ESEEM':None, 'ExpType':'5pDEER'}
        self.priorties = {'Auto': 150, 'MNR':300, 'Distance': 80, 'Single':200}

        self.priotityComboBox.addItems(list(self.priorties.keys()))
        self.correction_factor=1
        self.est_lambda = None

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
            elif model == 'ETH_AWG':
                from pyepr.hardware.ETH_awg import ETH_awg_interface
                self.spectromterInterface = ETH_awg_interface()
                self.spectromterInterface.savefolder = self.data_folder
                self.Bruker=False
                self.modeTuneDialog = ModeTune(self.spectromterInterface, gyro=self.gyro, threadpool=self.threadpool, current_folder=self.current_folder)
                self.modeTuneButton = QPushButton('Mode Tune')
                self.formLayout_2.addWidget(self.modeTuneButton)
                self.modeTuneButton.clicked.connect(self.modeTuneDialog.show)
                

            elif model == 'Bruker_MPFU':
                from pyepr.hardware.Bruker_MPFU import BrukerMPFU
                self.spectromterInterface = BrukerMPFU(filename_edit)
                self.spectromterInterface.savefolder = self.current_folder
                self.Bruker=True
            elif model == 'Bruker_AWG':
                from pyepr.hardware.Bruker_AWG import BrukerAWG
                self.spectromterInterface = BrukerAWG(filename_edit)
                self.spectromterInterface.savefolder = self.current_folder
                self.Bruker=True
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


    def select_resonator(self):
        key = self.resonatorComboBox.currentText()
        main_log.info(f"Selecting resonator {key}")
        self.LO = self.config['Resonators'][key]['Center Freq']
        self.fcDoubleSpinBox.setValue(self.LO)
        main_log.info(f"Setting LO to {self.LO} GHz")
    
    def change_LO(self):
        self.LO = self.fcDoubleSpinBox.value()
        main_log.info(f"Setting LO to {self.LO} GHz")

    @property
    def remaining_time(self):
        """
        Returns the remaining time in hours
        """
        return self.userinput['MaxTime']*3600 - (time.time() - self.starttime)

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

    def save_data(self,dataset,experiment,folder='data'):
        filename = create_save_name(self.userinput['sample'],experiment,True,self.userinput['project'],self.userinput['comment'])
        filename += ".h5"
        if folder == 'data':
            dataset.to_netcdf(os.path.join(self.data_folder,filename),engine='h5netcdf',invalid_netcdf=True)
        else:
            dataset.to_netcdf(os.path.join(self.current_folder,filename),engine='h5netcdf',invalid_netcdf=True)
        main_log.debug(f"Saving dataset to {filename}")



    def update_fieldsweep(self, dataset=None):

        if dataset is None:
            dataset = self.current_data['fieldsweep']
        else:
            self.current_data['fieldsweep'] = dataset
            self.save_data(dataset,'EDFS',folder='main')

        fsweep_analysis = epr.FieldSweepAnalysis(dataset)
        fsweep_analysis.calc_gyro()
        main_log.info(f"Calculated gyro {fsweep_analysis.gyro*1e3:.3f} MHz/T")
        if self.worker is not None:
            self.worker.set_noise_mode(np.min([fsweep_analysis.calc_noise_level(),20]))
        
        self.setup_plot_ax.cla()
        fsweep_analysis.plot(axs=self.setup_plot_ax,fig=self.setup_figure_canvas.figure)
        self.setup_figure_canvas.draw()
        self.Tab_widget.setCurrentIndex(1)

        # For Bruker recalculate d0
        if self.Bruker:
            main_log.info("Calculating d0 from Hahn Echo")
            B = self.LO/fsweep_analysis.gyro
            self.spectromterInterface.calc_d0_from_Hahn_Echo(B=B, LO=self.LO)

        if self.worker is not None:
            self.worker.update_gyro(fsweep_analysis.gyro)
        try:
            fit = self.config['autoDEER']['FieldSweep']['Fit']
        except KeyError:
            fit = False
        
        worker = Worker(fieldsweep_fit, fsweep_analysis,fit)
        worker.signals.result.connect(self.refresh_fieldsweep_after_fit)
        
        self.threadpool.start(worker)

 
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



    def update_respro(self, dataset=None):

        if dataset is None:
            dataset = self.current_data['respro']
        else:
            self.current_data['respro'] = dataset
            self.save_data(dataset,'ResPro',folder='main')
        
        f_lims = (self.config['Spectrometer']['Bridge']['Min Freq'], self.config['Spectrometer']['Bridge']['Max Freq'])

        # worker = Worker(respro_process, dataset, f_axis,self.current_results['fieldsweep'], cores=self.cores)
        worker = Worker(respro_process, dataset,f_lims, self.current_results['fieldsweep'], cores=self.cores)

        worker.signals.result.connect(self.refresh_respro)

        self.threadpool.start(worker)

    def create_setup_figure(self):
        fig, axs  = plt.subplots(1,1,figsize=(12.5, 6.28))
        self.setup_figure_canvas = FigureCanvas(fig)
        Navbar = NavigationToolbar2QT(self.setup_figure_canvas, self)
        Navbar.setMaximumHeight(24)
        self.respro_v_left.addWidget(self.setup_figure_canvas)
        self.respro_v_left.addWidget(Navbar)
        self.setup_plot_ax = axs


    def refresh_respro(self, *args):
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
            self.optimise_pulses_button()

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



    def optimise_pulses_button(self):
        if (self.pulses is None) or (self.pulses == {}):
            self.optimise_pulses()
        elif self.pulses['pump_pulse'] is None:
            self.optimise_pulses()
        else:
            self.update_optimise_pulses_figure()
        

    def optimise_pulses(self):
        resonator = self.current_results['respro']
        spectrum = self.current_results['fieldsweep']
        if self.pulses == {}:  # No pulses have been created yet
            # self.pulses = ad.build_default_pulses(self.AWG,tp = self.Min_tp)
            self.pulses = ad.create_pulses_shape(resonatorProfile=resonator,spectrum=spectrum)
        else:  # Reoptimise pulses
            if 'quickdeer' in self.current_results:
                r_min = 4.0
                max_tp = 256 # TODO
                pump_pulse_tp = self.pulses['pump_pulse'].tp.value
                if pump_pulse_tp > max_tp:
                    self.pulses = ad.create_pulses_shape(resonatorProfile=resonator,spectrum=spectrum,r_min=r_min)
                else:
                    main_log.info(f"Pulse are already optimised")
                    if self.waitCondition is not None: # Wake up the runner thread
                        self.waitCondition.wakeAll()
                    self.update_optimise_pulses_figure()
                    self.Tab_widget.setCurrentIndex(1)
        
        self.est_lambda = ad.calc_est_modulation_depth(spectrum, **self.pulses, respro=resonator)
        
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
        
        if self.waitCondition is not None: # Wake up the runner thread
            self.waitCondition.wakeAll()

        self.update_optimise_pulses_figure()
        self.Tab_widget.setCurrentIndex(1)
    
    def update_optimise_pulses_figure(self):
        pump_pulse = self.pulses['pump_pulse']
        exc_pulse = self.pulses['exc_pulse']
        ref_pulse = self.pulses['ref_pulse']

        self.setup_plot_ax.cla()
        ad.plot_overlap(self.current_results['fieldsweep'], pump_pulse, exc_pulse,ref_pulse, axs=self.setup_plot_ax,fig=self.setup_figure_canvas.figure, respro=self.current_results['respro'])
        self.setup_figure_canvas.draw()

        # update the pulse parameter grid
        type_to_pulse_hash = {epr.RectPulse:'Rect', epr.ChirpPulse:'Chirp', epr.HSPulse:'HS'}
        if isinstance(exc_pulse, epr.RectPulse):
            self.ExcFreqBox.setValue(param_in_MHz(exc_pulse.freq))
            self.ExcBWBox.setValue(exc_pulse.tp.value)
            self.ExcBWBox.setSuffix(' ns')
            self.ExcTypeLine.setText('Rect')
        else:
            center_freq = (param_in_MHz(exc_pulse.final_freq) + param_in_MHz(exc_pulse.init_freq))/2
            self.ExcFreqBox.setValue(center_freq)
            self.ExcBWBox.setValue(param_in_MHz(exc_pulse.bandwidth))
            self.ExcBWBox.setSuffix(' MHz')
            self.ExcTypeLine.setText(type_to_pulse_hash[type(exc_pulse)])
        
        if isinstance(ref_pulse, epr.RectPulse):
            self.RefFreqBox.setValue(param_in_MHz(ref_pulse.freq))
            self.RefBWBox.setValue(ref_pulse.tp.value)
            self.RefBWBox.setSuffix(' ns')
            self.RefTypeLine.setText('Rect')
        else:
            center_freq = (param_in_MHz(ref_pulse.final_freq) + param_in_MHz(ref_pulse.init_freq))/2
            self.RefFreqBox.setValue(center_freq)
            self.RefBWBox.setValue(param_in_MHz(ref_pulse.bandwidth))
            self.RefBWBox.setSuffix(' MHz')
            self.RefTypeLine.setText(type_to_pulse_hash[type(ref_pulse)])
        
        if isinstance(pump_pulse, epr.RectPulse):
            self.PumpFreqBox.setValue(param_in_MHz(pump_pulse.freq))
            self.PumpBWBox.setValue(pump_pulse.tp.value)
            self.PumpBWBox.setSuffix(' ns')
            self.PumpTypeLine.setText('Rect')
        else:
            center_freq = (param_in_MHz(pump_pulse.final_freq) + param_in_MHz(pump_pulse.init_freq))/2
            self.PumpFreqBox.setValue(center_freq)
            self.PumpBWBox.setValue(param_in_MHz(pump_pulse.bandwidth))
            self.PumpBWBox.setSuffix(' MHz')
            self.PumpTypeLine.setText(type_to_pulse_hash[type(pump_pulse)])

    def update_relax(self, dataset=None):
        if dataset is None:
            dataset = self.current_data['relax']
        else:
            dataset = dataset.epr.correctphasefull
            if 'relax' in self.current_data:
                # attempt to merge datasets
                main_log.info('Merging relax datasets')
                dataset = self.current_data['relax'].epr.merge(dataset)
            
            self.current_data['relax'] = dataset
            self.save_data(dataset, 'CP', folder='main')

        worker = Worker(relax_process, dataset)
        worker.signals.result.connect(self.refresh_relax)

        self.threadpool.start(worker)

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
            ax1 = fig.add_subplot(gs[0, 0])
            ax2 = fig.add_subplot(gs[1, 0])
            ax3 = fig.add_subplot(gs[1, 1])
            ax4 = fig.add_subplot(gs[0, 1])

            self.relax_ax = [ax1, ax2, ax3, ax4]

    def refresh_relax_figure(self):
    
        fig = self.relax_canvas.figure
        self.relax_ax[0].cla()
        relax1D_results = []
        if 'relax' in self.current_results:
            relax1D_results.append(self.current_results['relax'])
        if 'T2_relax' in self.current_results:
            relax1D_results.append(self.current_results['T2_relax'])

        epr.plot_1Drelax(*relax1D_results, axs=self.relax_ax[0], fig=fig, cmap=ad.primary_colors)
            
        if 'relax2D' in self.current_results:
            self.relax_ax[1].cla()
            self.current_results['relax2D'].plot2D(axs=self.relax_ax[3], fig=fig)

        self.relax_canvas.draw()


        
    def refresh_relax(self, fitresult):
        self.current_results['relax'] = fitresult


        self.refresh_relax_figure()
        self.label_eff = self.userinput['label_eff'] / 100

        if self.est_lambda is None:
            self.est_lambda = 0.4 
        self.aim_time = 2
        self.aim_MNR = 20
        
        self.initialise_deer_settings()

        # self.current_results['relax'].max_tau = max_tau
        # self.DipolarEvoMax.setValue(max_tau)
        # self.DipolarEvo2hrs.setValue(tau2hrs)
        self.Tab_widget.setCurrentIndex(2)
        
        if self.check_CP(fitresult): # CP decays passes test then it can proceed
            
            if self.worker is not None:
                # CP_decay = fitresult.func(fitresult.axis, *fitresult.fit_result[0]).data
                CP_decay = fitresult.fit_result.evaluate(fitresult.fit_model, fitresult.axis)*fitresult.fit_result.scale
                # Find the index when CP_decay is below 0.05
                CP_decay = CP_decay/CP_decay[0]
                CP_decay_bool = CP_decay < 0.05
                CP_decay_idx = np.where(CP_decay_bool)[0]
                if len(CP_decay_idx) == 0:
                    CP_decay_idx = len(CP_decay)
                else:
                    CP_decay_idx = CP_decay_idx[0]
                max_tau = fitresult.axis[CP_decay_idx]
                max_tau = epr.round_step(max_tau,1)
                main_log.info(f"Max tau {max_tau:.2f} us")
                self.worker.set_2D_max_tau(max_tau*2)

            if self.waitCondition is not None: # Wake up the runner thread
                self.waitCondition.wakeAll()

    def initialise_deer_settings(self):
        
        
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
        if 'relax' in self.current_results:
            relax_data['CP'] = self.current_results['relax']
        if 'T2_relax' in self.current_results:
            relax_data['Tm'] = self.current_results['T2_relax']
        if 'relax2D' in self.current_results:
            relax_data['Ref2D'] = self.current_results['relax2D']    
        
        #debug only, remove later
        store_pickle(relax_data,os.path.join(self.data_folder,'relax_data.pkl'))
        main_log.info(f'Calculating DEER settings with aim SNR {aim_SNR:.2f}, aim time {self.aim_time}hrs')
        self.deer_settings = ad.calc_DEER_settings(relax_data,exp, self.aim_time, aim_SNR,self.waveform_precision)
        
        self.deer_settings['criteria'] = self.aim_MNR

        self.deer_settings['autoStop'] = self.Time_autoStop_checkbox.isChecked()
        self.worker.update_deersettings(self.deer_settings)
        self.update_tau_delays_figure([aim_SNR],[self.aim_time],labels=[f"MNR = {self.aim_MNR}"])

        
        main_log.info(f"tau1 set to {self.deer_settings['tau1']:.2f} us")
        main_log.info(f"tau2 set to {self.deer_settings['tau2']:.2f} us")
        main_log.info(f"DEER Sequence set to {self.deer_settings['ExpType']}")

        return self.deer_settings

    def update_deer_settings(self,remaining_time=None):
        
        data = self.current_results['quickdeer']
        rec_tau = self.current_results['quickdeer'].rec_tau_max
        dt = self.current_results['quickdeer'].rec_dt * 1e3
        dt = epr.round_step(dt,self.waveform_precision)
        dt= 8
        mod_depth = data.MNR * data.noiselvl
        if remaining_time is None:
            remaining_time = self.remaining_time
            main_log.debug(f"Remaining time {remaining_time:.2f} hours")
        else:
            remaining_time = remaining_time
            main_log.debug(f"Measuring DEER for {remaining_time:.2f} hours")

        if self.deer_settings['ExpType'] == '4pDEER':
            relax = self.current_results['T2_relax']
        elif self.deer_settings['ExpType'] == '5pDEER':
            relax = self.current_results['relax']
        
        self.correction_factor = ad.calc_correction_factor(relax,self.current_results['quickdeer'])
        main_log.info(f"Correction factor {self.correction_factor:.3f}")

        # Assemble all relaxation data
        if (self.Exp_types.currentText() == '4pDEER'):
            exp = '4pDEER'
        else:
            exp = 'auto'

        if self.userinput['priority'].lower() == 'single':
            SNR_target = self.priorties[self.userinput['priority']]
            MNR_target = SNR_target
            single_mode = True
            exp = '4pDEER'
        else:
            MNR_target = self.priorties[self.userinput['priority']]
            SNR_target = MNR_target/(mod_depth)
        main_log.info(f"SNR target {SNR_target:.2f}")


        relax_data = {}
        if 'relax' in self.current_results:
            relax_data['CP'] = self.current_results['relax']
        if 'T2_relax' in self.current_results:
            relax_data['Tm'] = self.current_results['T2_relax']
        if 'relax2D' in self.current_results:
            relax_data['Ref2D'] = self.current_results['relax2D']         

        #debug only, remove later
        store_pickle(relax_data,os.path.join(self.data_folder,'relax_data.pkl'))

        # Calculate the optimal DEER settings using relaxation data

        self.deer_settings = ad.calc_DEER_settings(relax_data,exp,remaining_time,SNR_target,self.waveform_precision,corr_factor=self.correction_factor,rec_tau=rec_tau)

        # self.deer_settings['dt'] = dt
        self.deer_settings['criteria'] = MNR_target
        
        self.deer_settings['autoStop'] = self.Time_autoStop_checkbox.isChecked()
        self.worker.update_deersettings(self.deer_settings)
        self.update_tau_delays_figure([SNR_target],[remaining_time],labels=[f"MNR = {MNR_target}"])
        
        main_log.info(f"tau1 set to {self.deer_settings['tau1']:.2f} us")
        main_log.info(f"tau2 set to {self.deer_settings['tau2']:.2f} us")
        main_log.info(f"DEER Sequence set to {self.deer_settings['ExpType']}")

    def update_tau_delays_figure(self, SNRs, MeasTimes, labels=None):

        fig = self.relax_canvas.figure
        axs = self.relax_ax[2]
        axs.cla()
        
        if 'relax' in self.current_results:
            CP_analysis = self.current_results['relax']
            ad.plot_optimal_tau(CP_analysis,SNRs,MeasTimes,MaxMeasTime=36, labels=['5pDEER'],fig=fig,axs=axs,cmap=[epr.primary_colors[0]],corr_factor=self.correction_factor)

        if 'relax2D' in self.current_results:
            ad.plot_optimal_tau(self.current_results['relax2D'],SNRs,MeasTimes,MaxMeasTime=36, labels=['4pDEER'],fig=fig,axs=axs,cmap=[epr.primary_colors[1]],corr_factor=self.correction_factor)
        elif 'T2_relax' in self.current_results:
            ad.plot_optimal_tau(self.current_results['T2_relax'],SNRs,MeasTimes,MaxMeasTime=36, labels=['4pDEER'],fig=fig,axs=axs,cmap=[epr.primary_colors[1]],corr_factor=self.correction_factor)

        axs.set_title(labels[0])

    def update_relax2D(self, dataset=None):
        if dataset is None:
            dataset = self.current_data['relax2D']
        else:
            self.current_data['relax2D'] = dataset
            self.save_data(dataset,'2D_DEC',folder='main')

        # Since there is no fitting in the 2D data analysis it can be run in the main threads

        relax2DAnalysis = ad.RefocusedEcho2DAnalysis(dataset)
        self.current_results['relax2D'] = relax2DAnalysis

        self.refresh_relax_figure()

        self.initialise_deer_settings()

        if self.waitCondition is not None: # Wake up the runner thread
            self.waitCondition.wakeAll()
        

    
    def update_T2(self,dataset=None):
        if dataset is None:
            dataset = self.current_data['T2_relax']
        else:
            dataset = dataset.epr.correctphasefull
            if 'T2_relax' in self.current_data:
                # attempt to merge datasets
                main_log.info('Merging relax datasets')
                dataset = self.current_data['T2_relax'].epr.merge(dataset)
            
            self.current_data['T2_relax'] = dataset
            self.save_data(dataset,'T2',folder='main')

        
        d_ESEEM = epr.detect_ESEEM(dataset,'deuteron')
        p_ESEEM = epr.detect_ESEEM(dataset,'proton')
        d_ESEEM = False # Turn off ESEEM detection
        p_ESEEM = False
        if d_ESEEM:
            self.deer_settings['ESEEM'] = 'deuteron'
            main_log.info(f"Detected deuteron ESEEM")
        elif p_ESEEM:
            self.deer_settings['ESEEM'] = 'proton'
            main_log.info(f"Detected proton ESEEM")
        else:
            self.deer_settings['ESEEM'] = None
            main_log.info(f"No ESEEM detected")
        # Since the T2 values are not used for anything there is no need pausing the spectrometer    

        T2_worker = Worker(T2_process, dataset)
        T2_worker.signals.result.connect(self.refresh_T2)
        # T2_worker.signals.result.connect(self.check_T2)

        self.threadpool.start(T2_worker)

    def check_T2(self, fitresult):
        # Check if the T2 measurment is too short. 

        test_result = fitresult.check_decay()
        test_dt = fitresult.axis[1].values - fitresult.axis[0].values
        test_dt *= 1e3

        if test_result == 0:
            # if self.waitCondition is not None: # Wake up the runner thread
            #     self.waitCondition.wakeAll()
            return True
        
        elif test_result == -1:  # The trace needs to be longer
            new_dt = epr.round_step(test_dt*2, self.waveform_precision)
        elif test_result == 1:  # The trace needs to be shorter
            new_dt = epr.round_step(test_dt/2, self.waveform_precision)

        new_tmin = fitresult.axis[-1].values
        new_tmin += new_dt*1e-3
        nAvgs = fitresult.dataset.attrs['nAvgs']

        self.initialise_deer_settings()
        
        if self.worker is not None:
            self.worker.run_T2_relax(dt=new_dt, tmin=new_tmin, averages=nAvgs, autoStop=False)
            return False
        else:
            return True

    def check_CP(self, fitresult):
        # Check if the CP measurment is too short. 

        test_result = fitresult.check_decay()
        main_log.debug(f"CP test result {test_result}")

        test_dt = fitresult.axis[1].values - fitresult.axis[0].values
        test_dt *= 1e3
        if test_result == 0:
            # if self.waitCondition is not None: # Wake up the runner thread
            #     self.waitCondition.wakeAll()
            return True
        elif test_result == -1: # The trace needs to be longer
            new_dt = epr.round_step(test_dt*2, self.waveform_precision)
        elif test_result == 1: # The trace needs to be shorter
            new_dt = epr.round_step(test_dt/2, self.waveform_precision)

        new_tmin = fitresult.axis[-1].values
        new_tmin += new_dt*1e-3
        nAvgs = fitresult.dataset.attrs['nAvgs']
        
        if self.worker is not None:
            self.worker.run_CP_relax(dt=new_dt, tmin=new_tmin*2, averages=nAvgs, autoStop=False)
            return False
        else:
            return True

    def refresh_T2(self, fitresult):
        self.current_results['T2_relax'] = fitresult
        self.refresh_relax_figure()

        if self.check_T2(fitresult):

            if self.waitCondition is not None:  # Wake up the runner thread
                self.waitCondition.wakeAll()

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
                time = np.min([self.aim_time,self.remaining_time])
                self.worker.repeat_quickdeer()
            else:
                time = None

            self.update_deer_settings(remaining_time=time)

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

    def update_reptime(self, dataset=None):
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
            ReptimeRecovery = 0.75
        opt_reptime = reptime_analysis.calc_optimal_reptime(ReptimeRecovery)

        if (opt_reptime*1e-3 > 8) or (opt_reptime*1e-3 < 0.5):
            main_log.warning(f"Reptime optimisation failed. Setting to default for spin system")
            opt_reptime = 3e3

        self.current_results['reptime'] = reptime_analysis
        if self.worker is not None:
            self.worker.update_reptime(opt_reptime)
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
        userinput['SampleConc'] = SampleConcComboBox_opts[self.SampleConcComboBox.currentText()]
        userinput['tp'] = self.Min_tp

        self.userinput = userinput

        if self.priotityComboBox.currentText().lower() == 'single':
            self.operating_mode = 'single'
        else:
            self.operating_mode = self.Exp_types.currentText()

        if advanced:
                        
            if self.Tau1Value.value() > 0 or self.Tau2Value.value() > 0:
                self.fixed_tau = {'tau1':self.Tau1Value.value(),'tau2':self.Tau2Value.value(),'tau3':self.Tau3Value.value()}
                # Skip the relaxation data analysis and go straight to DEER
                if self.Exp_types.currentText() == 'auto':
                    self.deer_settings['ExpType'] = '5pDEER'
                elif self.Exp_types.currentText() == 'Ref2D':
                    print("ERROR: Ref2D not supported with fixed tau")
                    return None
                else:
                    self.deer_settings['ExpType'] = self.Exp_types.currentText()
                
                self.deer_settings['tau1'] = self.Tau1Value.value()
                self.deer_settings['tau2'] = self.Tau2Value.value()
                self.deer_settings['tau3'] = self.Tau3Value.value()
                self.deer_settings['criteria'] = self.priorties[self.userinput['priority']]
                self.deer_settings['autoStop'] = self.Time_autoStop_checkbox.isChecked()
                self.deer_settings = calc_dt_from_tau(self.deer_settings)
                main_log.info(f"tau1 set to {self.deer_settings['tau1']:.2f} us")
                main_log.info(f"tau2 set to {self.deer_settings['tau2']:.2f} us")
                main_log.info(f"DEER Sequence set to {self.deer_settings['ExpType']}")
            
            else:
                self.fixed_tau = None
            

        
        
        if not 'ESEEM' in self.deer_settings:
            self.deer_settings['ESEEM'] = None

        if self.Exp_types.currentText() in ['4pDEER','Ref2D']:
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
        
        self.starttime = time.time()

        self.worker.update_deersettings(self.deer_settings)
    
        self.worker.signals.status.connect(self.msgbar.setText)
        self.worker.signals.status.connect(main_log.info)
        self.worker.signals.fsweep_result.connect(self.update_fieldsweep)
        self.worker.signals.fsweep_result.connect(lambda x: self.save_data(x,'EDFS'))
        self.worker.signals.respro_result.connect(self.update_respro)
        self.worker.signals.respro_result.connect(lambda x: self.save_data(x,'ResPro'))
        # self.worker.signals.optimise_pulses.connect(self.optimise_pulses)
        self.worker.signals.relax_result.connect(self.update_relax)
        self.worker.signals.relax_result.connect(lambda x: self.save_data(x,'CP'))
        self.worker.signals.T2_result.connect(self.update_T2)
        self.worker.signals.T2_result.connect(lambda x: self.save_data(x,'T2'))

        self.worker.signals.Relax2D_result.connect(self.update_relax2D)
        self.worker.signals.Relax2D_result.connect(lambda x: self.save_data(x,'2D_DEC'))

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
    
    def RunFullyAutoDEER(self):
        return self.RunAutoDEER(advanced=False)

    def RunAdvancedAutoDEER(self):
        return self.RunAutoDEER(advanced=True)

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
    
