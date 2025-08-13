from PyQt6.QtWidgets import QDialog
from PyQt6 import uic
import PyQt6.QtCore as QtCore

from pyepr import HahnEchoSequence, RectPulse, Parameter, Interface
import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvas, NavigationToolbar2QT
import matplotlib.pyplot as plt
import time
import os
from autodeer.gui.tools import WorkerSignals
from autodeer.colors import primary_colors

import datetime
QtCore.QDir.addSearchPath('icons', 'gui/resources')
package_directory = os.path.dirname(os.path.abspath(__file__))

from scipy.optimize import curve_fit


def save_data(dataset,filename,folder=''):
        timestamp = datetime.datetime.now().strftime(r'%Y%m%d_%H%M_')
        fullname = timestamp + '_' +filename +'.h5'
        dataset.to_netcdf(os.path.join(folder,fullname),engine='h5netcdf',invalid_netcdf=True)

def phasecorrect_all_points(dataset):
    dta_angle = np.angle(dataset.data)
    data = dataset.data * np.exp(-1j * dta_angle)
    new_dataset = dataset.copy()
    new_dataset.data = data
    return new_dataset

class ModeTune(QDialog):
    dataUpdated = QtCore.pyqtSignal(dict)

    def __init__(self, interface: Interface, gyro:float = 0.002803632236095, threadpool= None,current_folder=''):
        super().__init__()
        self.interface = interface
        self.ui = uic.loadUi(os.path.join(package_directory, 'modetuneUI.ui'), self)

        self.create_figure()
        self.ui.gyro_ratio.setValue(gyro *1e3)
        if threadpool is not None:
            self.threadpool = threadpool
        else:
            self.threadpool = QtCore.QThreadPool()

        if current_folder  == '':
            self.current_folder = os.getcwd()
        else:
            self.current_folder = current_folder

        self.ui.searchbutton.clicked.connect(self.start_mode_tune)

    def start_mode_tune(self):
        self.mode_ax.cla()
        gyro_ratio = self.ui.gyro_ratio.value()/1e3
        center_freq = self.ui.center_freq.value()
        scan_range = self.ui.scan_range.value() / 1e3

        n_points = 50
        scan_step  = scan_range/n_points

        freq = Parameter("freq", -scan_range/2, step=scan_step,dim=n_points, unit="GHz")
        pi2_pulse = RectPulse(tp=16,freq=freq, scale=0.1,flipangle=np.pi/2)
        pi_pulse = RectPulse(tp=32,freq=freq, scale=0.1,flipangle=np.pi)
        B = Parameter(
            "B",((center_freq)/gyro_ratio), start=-(scan_range/2)/gyro_ratio, step=scan_step/gyro_ratio, dim=n_points,
            unit="Guass",link=freq,description="B0 Field" )
        
        seq = HahnEchoSequence(freq=center_freq, B = B, tau=500, reptime=3e3, shots=100, averages = 1, pi2_pulse=pi2_pulse, pi_pulse=pi_pulse)
        seq.pulses[0].freq = freq
        seq.pulses[1].freq = freq
        seq.pulses[2].freq = freq
        seq.pulses[2].tp.value = 512
        seq.evolution([freq])

        self.interface._launch(seq,savename="modetune",IFgain=0,mode='tune')

        # Start thread
        self.worker = get_dataWorker(self.interface)
        self.worker.signals.result.connect(self.update_figure)
        self.worker.signals.finished.connect(self.update_dip)

        self.threadpool.start(self.worker)


    def create_figure(self):
        fig, axs  = plt.subplots(1,1,figsize=(12.5, 6.28),layout='constrained')
        self.modetTunecanvas = FigureCanvas(fig)
        self.figure_v_layout.addWidget(self.modetTunecanvas)
        Navbar = NavigationToolbar2QT(self.modetTunecanvas, self)
        Navbar.setMaximumHeight(20)
        self.figure_v_layout.addWidget(Navbar)
        self.mode_ax = axs
        self.mode_ax.set_ylim(0,1)
        self.mode_ax.set_xlabel("Frequency (GHz)")
        self.mode_ax.set_ylabel("Intensity (a.u.)")
        self.modetTunecanvas.draw()


    def update_figure(self,dataset,done=False):
        self.dataset = dataset
        fig = self.modetTunecanvas.figure
        self.mode_ax.cla()

        dataset_pc = phasecorrect_all_points(dataset)
        dataset_pc /= np.max(dataset_pc.data)
        self.mode_ax.set_ylim(0,1)
        self.mode_ax.plot(dataset_pc.freq + dataset_pc.pulse0_freq, dataset_pc.data)
        self.mode_ax.set_xlabel("Frequency (GHz)")
        self.mode_ax.set_ylabel("Intensity (a.u.)")
        save_data(dataset,"modetune",folder=self.current_folder)
        self.modetTunecanvas.draw()

    def update_dip(self):
        dataset = self.dataset
        dataset_pc = phasecorrect_all_points(dataset)
        dataset_pc /= np.max(dataset_pc.data)

        # fit with gaussian
        x = dataset_pc.freq + dataset_pc.pulse0_freq
        y = dataset_pc.data

        self.mode_ax.cla()
        self.mode_ax.plot(x, y, '.', color='0.7', label='Data', ms=6)
        self.mode_ax.set_xlabel("Frequency (GHz)")
        self.mode_ax.set_ylabel("Intensity (a.u.)")

        def gaussian(x, A, mu, sigma):
            return A * np.exp(-((x - mu) ** 2) / (2 * sigma ** 2))
        try:
            popt, pcov = curve_fit(gaussian, x, y, p0=[1, dataset_pc.freq, 0.1])
        except:
            popt = None
       

        if popt is not None:
            fit_data = gaussian(x, *popt)
            # check r^2
            residuals = y - fit_data
            ss_res = np.sum(residuals ** 2)
            ss_tot = np.sum((y - np.mean(y)) ** 2)
            R2 = 1 - (ss_res / ss_tot)
            fc = popt[1]
            sigma = popt[2]
            Q = fc / sigma
            self.fc = fc
            self.Q = Q

            if R2 < 0.5:
                #set QDoubleSpinBox color to red
                self.ui.calc_freq.setStyleSheet("background-color: red")
            elif R2 < 0.75:
                #set QDoubleSpinBox color to yellow
                self.ui.calc_freq.setStyleSheet("background-color: yellow")
            else:
                #set QDoubleSpinBox color to green
                self.ui.calc_freq.setStyleSheet("background-color: green")

            self.ui.calc_freq.setValue(popt[1])
            self.ui.qDoubleSpinBox.setValue(Q)
            self.mode_ax.plot(x, gaussian(x, *popt), label='fit',color=primary_colors[0])

        # update plot
        self.mode_ax.legend()
        self.modetTunecanvas.draw()
        pass
   
    def send_data_to_main(self):
        # Prepare the data to send
        data = {
            "fc": self.fc,
            "Q": self.Q,
        }
        # Emit the signal with the data
        self.dataUpdated.emit(data)

class get_dataWorker(QtCore.QRunnable):

    result_signal = QtCore.pyqtSignal(object)
    finished_signal = QtCore.pyqtSignal()

    def __init__(self, interface: Interface):
        super(get_dataWorker, self).__init__()
        self.interface = interface
        self.signals = WorkerSignals()
    
    
    @QtCore.pyqtSlot()
    def run(self):
        while self.interface.isrunning():
            try:
                dataset = self.interface.acquire_dataset()
            except:
                continue
            self.signals.result.emit(dataset)
            time.sleep(2)
        dataset = self.interface.acquire_dataset()
        self.signals.result.emit(dataset)
        self.signals.finished.emit()

