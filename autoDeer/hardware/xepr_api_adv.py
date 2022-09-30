import numpy as np
import time
import os
import XeprAPI
from openepr import dataset

import logging

hw_log = logging.getLogger('hardware.Xepr')

# =============================================================================


# class dataset:
#     def __init__(self, time: np.ndarray, data: np.ndarray, 
#                  cur_exp: XeprAPI.Xepr.XeprExperiment = None) -> None:

#         self.time = time
#         self.data = data
#         self.dim: int = int(len(time))
#         self.size = np.shape(data)
#         if cur_exp is not None:
#             self.scans_done = cur_exp.getParam("NbScansDone").value
#             self.scans_todo = cur_exp.getParam("NbScansToDo").value
#             self.shrt = cur_exp.getParam("ShotRepTime").value
#             self.shot_p_point = cur_exp.getParam("ShotsPLoop").value
#         pass


# ============================================================================


hardware_meta = {  # This dictionary should be moved into a config file 
    "Type":             "Complete Spectrometer",
    "Manufacturer":     "Bruker",
    "Model":            "E600",
    "Local name":       "C_Floor",
    "Acq. Resolution":  2,
    "Pulse Resolution": 2,
    "AWG":              False,
    "Min Freq":         33,
    "Max Freq":         35,
 } 


class xepr_api:
    def __init__(self) -> None:
        self.Xepr = None
        self.cur_exp = None
        self.hidden = None
        self._tmp_dir = None
        self.XeprCmds = None
        self.spec_config = hardware_meta
        pass

    def connect(self) -> None:
        """
        Connect to the Xepr Spectrometer.
        """
        self.find_Xepr()
        self.find_cur_exp()
        self.find_hidden()
        pass

    def _set_Xepr_global(self, Xepr_inst):
        self.Xepr = Xepr_inst
        self.XeprCmds = self.Xepr.XeprCmds

    def _get_Xepr_global(self):
        if self.Xepr is not None:
            return self.Xepr
        else:
            raise RuntimeError("Can't find XEPR instance")

    def find_Xepr(self):
        open_xepr_instances = XeprAPI.getXeprInstances()

        if len(open_xepr_instances) < 1:
            raise RuntimeError("No Open Xepr API Instances, please open one by"
                               ":\n \"processing\" -> \"XeprAPI\" -> \"Enable "
                               "Xepr API\"")

        if self.Xepr is None:
            try:
                self.Xepr_local = XeprAPI.Xepr()
            except OSError:
                print('Xepr API: Could not connect to any Xepr instance.')
            except RuntimeError:
                raise RuntimeError("There is already a connection from Xepr to"
                                   "this python Kernel.\n Please use the "
                                   "correct python object or restart the "
                                   "kernel ")
            except:
                RuntimeError("Can't connect to Xepr: Please check Xepr is "
                             "running and open to API")
            else:
                self._set_Xepr_global(self.Xepr_local)
        else:
            print('Already Connected to Xepr!')

    def _set_cur_exp_global(self, cur_exp_inst):
        self.cur_exp = cur_exp_inst

    def _get_cur_exp_global(self):
        if self.cur_exp is not None:
            return self.cur_exp
        else:
            raise RuntimeError("Can't find current experiment")

    def find_cur_exp(self):
        """
        Try and finds the current experiment
        """
        if self.Xepr is None:
            raise RuntimeError("Please connect API to Xepr")

        try:
            self.cur_exp = self.Xepr.XeprExperiment()
        except:
            print("Can't find the current experiment. Attempting to load it")
            self.Xepr.XeprCmds.aqExpSelect("Experiment")
            try:
                self.cur_exp = self.Xepr.XeprExperiment()
            except:
                RuntimeError("Can't find the current experiment. Please create" 
                             "an experiment with name 'Experiment'")
            else:
                print("Experiment found")
        self._set_cur_exp_global(self.cur_exp)
        return self.cur_exp

    def find_hidden(self):
        if self.Xepr is not None:
            self.hidden = self.Xepr.XeprExperiment("AcqHidden")

    def is_exp_running(self):
        return self.cur_exp.isRunning

    def acquire_dataset(self) -> dataset:
        """
        This function acquire the dataset, this work both for a running 
        experiment or once it has finished.
        """
        dataclass = self.Xepr.XeprDataset()
        size = dataclass.size
        data_dim = len(size)
        data = dataclass.O
        params = {
            "scans_done": self.cur_exp.getParam("NbScansDone").value,
            "scans_todo": self.cur_exp.getParam("NbScansToDo").value,
            "shrt": self.cur_exp.getParam("ShotRepTime").value,
            "shots": self.cur_exp.getParam("ShotsPLoop").value
            }

        if data_dim == 1:
            # We have a 1D dataset
            t = dataclass.X
            hw_log.debug('Acquired Dataset')
            return dataset(t, data, params)
        elif data_dim == 2:
            # we have a 2D dataset
            t1 = dataclass.X
            t2 = dataclass.Y
            hw_log.debug('Acquired Dataset')
            return dataset([t1, t2], data, params)

    def acquire_scan(self):
        """
        This script detects the end of the scan and acquires the data set. 
        This requires that the experiment is still running, or recently 
        finished. Once it has been saved this will no longer work.
        """
        if self.is_exp_running():
        
            current_scan = self.cur_exp.getParam("NbScansDone").value
            # x_length = int(self.cur_exp.getParam("XSpecRes").value)
            time_per_point = self.cur_exp.getParam("ShotRepTime").value * 1e-6\
                * self.cur_exp.getParam("ShotsPLoop").value*2
            # trace = np.zeros(x_length, dtype=np.complex64)
            while self.cur_exp.getParam("NbScansDone").value == current_scan:
                time.sleep(time_per_point)
            time.sleep(time_per_point)
            return self.acquire_dataset()
        else:
            return self.acquire_dataset()

    def acquire_scan_at(self, scan_num: int):
        """
        This script acquires the scan after a specific number of scans
        """
        x_length = int(self.cur_exp.getParam("XSpecRes").value)
        time_per_point = self.cur_exp.getParam("ShotRepTime").value * 1e-6 \
            * self.cur_exp.getParam("ShotsPLoop").value * 2
        while self.cur_exp.getParam("NbScansDone").value != scan_num:
            time.sleep(time_per_point * x_length / 2)
        return self.acquire_scan()

    def acquire_scan_2d(self):
        """
        This function acquires the dataset after a full 2D scan.
        This is done by identifying the number of scan steps per sweep and 
        acquiring the data on that scan.
        This requires that the experiment has not been saved. 
        """
        if self.is_exp_running():
            total_num_scan = self.cur_exp.getParam("NbScansToDo").value
            total_num_sweeps = self.cur_exp.getParam("SweepsPExp").value
            scans_per_sweep = total_num_scan/total_num_sweeps

            if not scans_per_sweep.is_integer():
                raise RuntimeError("There is a non integer number of scans per"
                                   " sweep")

            current_scan = self.cur_exp.getParam("NbScansDone").value
            current_sweep = np.floor(current_scan/scans_per_sweep)
            next_scan_target = (current_sweep + 1) * scans_per_sweep

            return self.acquire_scan_at(next_scan_target)
        else:
            return self.acquire_scan()

    def set_PulseSpel_var(self, variable: str, value: int):
        """
        This can be used to change any pulse spell variable whilst the 
        experiment is running. These changes take effect at the beginning of 
        the next scan
        """
        hw_log.debug(f'Set pulse spell var {str(variable)} to {int(value)}')
        self.cur_exp["ftEpr.PlsSPELSetVar"].value = str(variable) + " = " + \
            str(int(value))

    def set_ReplaceMode(self, state=False):
        if state:
            value = 'On'
            hw_log.warning('Replace mode turned on')
            print('DANGER: Replace Mode turned ON!')
        else:
            value = 'Off'
            hw_log.info('Replace mode turned off')

        self.cur_exp["ftEpr.ReplaceMode"].value = value

    def set_PhaseCycle(self, state=True):
        if state:
            hw_log.info("On-Board Phase Cycling turned on")
            self.cur_exp["PCycleOn"].value = state
        else:
            hw_log.info("On-Board Phase Cycling turned off")
            self.cur_exp["PCycleOn"].value = state
  
        return state

    def get_PulseSpel_exp_filename(self):
        return os.path.basename(self.cur_exp["ftEPR.PlsSPELPrgPaF"].value)

    def get_PulseSpel_exp_filepath(self):
        self.save_PulseSpel_exp()
        return self.cur_exp["ftEPR.PlsSP§ELPrgPaF"].value

    def set_PulseSpel_exp_filepath(self, fullpath):
        self.save_PulseSpel_exp()
        self.Xepr.XeprCmds.aqPgLoad(fullpath)

    def get_PulseSpel_def_filename(self):
        return os.path.basename(self.cur_exp["ftEPR.PlsSPELGlbPaF"].value) 

    def get_PulseSpel_def_filepath(self):
        return self.cur_exp["ftEPR.PlsSPELGlbPaF"].value

    def set_PulseSpel_def_filepath(self, fullpath):
        self.save_PulseSpel_def()
        self.Xepr.XeprCmds.aqPgDefLoad(fullpath)

    def get_PulseSpel_phase_cycling(self):
        return self.cur_exp["ftEPR.PlsSPELLISTSlct"].value

    def set_PulseSpel_phase_cycling(self, name):
        self.cur_exp["ftEPR.PlsSPELLISTSlct"].value = name

        if self.cur_exp["ftEPR.PlsSPELLISTSlct"].value != name:
            print("WARNING: Phase cycling did not change")
            return 0
        else:
            return 1

    def get_PulseSpel_experiment(self):
        return self.cur_exp["ftEPR.PlsSPELEXPSlct"].value

    def set_PulseSpel_experiment(self, name):
        self.cur_exp["ftEPR.PlsSPELEXPSlct"].value = name

        if self.cur_exp["ftEPR.PlsSPELEXPSlct"].value != name:
            print("WARNING: Pulse Spel Experiment did not change")
            return 0
        else:
            return 1

    def save_PulseSpel_exp(self, name: str = None):
        if name is None:
            # Save as a temp file
            if self._tmp_dir is None:
                try:
                    os.mkdir('/tmp/autoDeer/')
                    self._tmp_dir = '/tmp/autoDeer/'
                except:
                    self._tmp_dir = '/tmp/autoDeer/'
                    print("temp directory already exists")
            timestr = time.strftime("%Y%m%d-%H%M%S")
            path = self._tmp_dir + "pulsespel_" + timestr + ".exp"
        else:
            path = name
        hw_log.debug(f'Saved Pulse Spel experiment to: {path}')
        self.Xepr.XeprCmds.aqPgSaveAs(path)

    def save_PulseSpel_def(self, name=None):
        if name is None:
            # Save as a temp file
            if self._tmp_dir is None:
                try:
                    os.mkdir('/tmp/autoDeer/')
                    self._tmp_dir = '/tmp/autoDeer/'
                except:
                    self._tmp_dir = '/tmp/autoDeer/'
                    print("temp directory already exists")
            timestr = time.strftime("%Y%m%d-%H%M%S")
            path = self._tmp_dir + "pulsespel_" + timestr + ".def"
        else:
            path = name
        hw_log.debug(f'Saved Pulse Spel definition to: {path}')
        self.Xepr.XeprCmds.aqPgDefSaveAs(path)

    def get_Acquisition_mode(self):
        return self.cur_exp["ftEPR.FTAcqModeSlct"].value

    def set_Acquisition_mode(self, mode: int):
        """mode=0: Run from tables, mode=1: Run from Pulse Spel,
        mode=2:Read transient, mode=3:Start Transient"""
        if mode == 0:
            self.cur_exp["ftEPR.FTAcqModeSlct"].value = 'Run from Tables'
        elif mode == 1:
            self.cur_exp["ftEPR.FTAcqModeSlct"].value = 'Run from PulseSPEL'
        elif mode == 2:
            self.cur_exp["ftEPR.FTAcqModeSlct"].value = 'Read Transient'
        elif mode == 3:
            self.cur_exp["ftEPR.FTAcqModeSlct"].value = 'Start Transient'
        else:
            print('Acquisition Mode not changed. Input error.')
            return 0

    def compile_PulseSpel_prg(self):
        self.Xepr.XeprCmds.aqPgShowPrg()
        self.Xepr.XeprCmds.aqPgCompile()
        time.sleep(0.5)
        pass

    def compile_PulseSpel_def(self):
        self.Xepr.XeprCmds.aqPgShowDef()
        self.Xepr.XeprCmds.aqPgCompile()
        time.sleep(0.5)
        pass

    def run_exp(self):
        self.cur_exp.aqExpRun()
        hw_log.info('Experiment started')
        time.sleep(5)  # Required
        pass

    def stop_exp(self):
        self.cur_exp.aqExpStop()
        hw_log.info('Experiment stopped')
        pass

    def abort_exp(self):
        self.cur_exp.aqExpAbort()
        hw_log.info('Experiment aborted')
        pass

    def xepr_save(self, path: str, title: str = None):
        """_summary_

        Parameters
        ----------
        path : str
            _description_
        title : str, optional
            _description_, by default None
        """
        xepr_file_limit = 70
        directory, basename = os.path.split(path)
        if not title:
            title = basename
        if not os.path.exists(directory):
            os.makedirs(directory)

        if len(basename) > xepr_file_limit:
            print("File name is too long. The file_name will be truncated, but"
                  " the full name will be saved as the title")
            path = directory + "/" + basename[:70]

        self.Xepr.XeprCmds.vpSave("Current Primary", title, path)
        hw_log.info(f'Saved data to: {path}')
        time.sleep(0.5)
        self.Xepr.XeprCmds.aqExpSelect("Experiment")
        time.sleep(0.5)
        pass

    def dataset_save(self):
        pass

    def get_field(self) -> int:
        """ This returns the central field"""
        return self.cur_exp['CenterField'].value

    def set_field(self, val: int, hold: bool = True) -> int:
        """ This sets the central field"""
        self.cur_exp['CenterField'].value = val
        time.sleep(2)  # Always wait 2s after a field change
        hw_log.info(f'Field position set to {val} G')
        if hold is True:
            while self.cur_exp['FieldWait'] is True:
                time.sleep(0.5)
        return self.cur_exp['CenterField'].value

    def get_counterfreq(self) -> float:
        """ This returns the current freq counter"""
        return self.cur_exp['FrequencyMon'].value

    def set_sweep_width(self, val: int) -> int:
        self.cur_exp['SweepWidth'].value = val
        hw_log.info('Field sweep width set to {val} G')
        return self.cur_exp['SweepWidth'].value
    
    def get_sweep_width(self) -> int:
        return self.cur_exp['SweepWidth'].value

    def set_freq(self, val: np.float128) -> float:
        """
        This sets bridge frequency, and works through a polynomial
        approximation. This might need to be adjusted
        for different spectrometers
        """
        f_pol = [1.420009750632201e4,
                 -5.118516287710228e3,
                 2.092103562165744e2,
                 -2.034307248428457,
                 0, 0]

        pol_func = np.polynomial.polynomial.Polynomial(f_pol)
        pos = round(pol_func(val))
        self.hidden['Frequency'].value = pos
        hw_log.info(f'Frequency set to {val} at position {pos}')
        return self.hidden['Frequency'].value

    def get_freq(self) -> float:
        """ This returns the current bridge frequency"""
        return self.hidden['Frequency'].value

    def get_spec_config(self) -> str:
        """get_spec_config Gets the name of the current spectrometer
        configuration file.

        Returns
        -------
        str
            Returns the name of the current spectrometer configuration file
        """

        return self.hidden['PlsPrgCalDbName'].value

    def set_spec_config(self, name: str = 'Normal') -> str:
        """set_spec_config Sets the name of the current spectrometer
        configuration file

        Parameters
        ----------
        name : str, optional
            The file name of config file. Normal and AWG, shortcut to the
            standard types, by default 'Normal'

        Returns
        -------
        str
            Returns the name of the current spectrometer configuration file
        """
        if name == 'Normal':
            config_file = 'Q_TWT_Jun21'
        elif name == 'AWG':
            config_file = 'Q_awgins2013'
        else:
            config_file = name
        
        if self.get_spec_config() != config_file:
            self.hidden['PlsPrgCalDbName'].value = config_file
            self.hidden['PlsPrgCalDbLoad']
            self.hidden['ApplyCfg']

    def get_attenuator(self, channel: str) -> float:
        if channel == 'Main':
            atten_channel = 'PowerAtten'
        elif channel == '+<x>':
            atten_channel = 'BrXAmp'
        elif channel == '-<x>':
            atten_channel = 'BrMinXAmp'
        elif channel == '+<y>':
            atten_channel = 'BrYAmp'
        elif channel == '-<y>':
            atten_channel = 'BrMinYAmp'
        elif channel == 'ELDOR':
            atten_channel = 'ELDORAtt'

        return self.api.hidden[atten_channel].value

    def set_attenuator(self, channel: str, value) -> float:
        if channel == 'Main':
            atten_channel = 'PowerAtten'
        elif channel == '+<x>':
            atten_channel = 'BrXAmp'
        elif channel == '-<x>':
            atten_channel = 'BrMinXAmp'
        elif channel == '+<y>':
            atten_channel = 'BrYAmp'
        elif channel == '-<y>':
            atten_channel = 'BrMinYAmp'
        elif channel == 'ELDOR':
            atten_channel = 'ELDORAtt'
        
        self.api.hidden[atten_channel].value = value

        return self.get_attenuator(channel)

    def get_phase(self, channel: str) -> float:
        if channel == 'Main':
            phase_channel = 'PowerAtten'
        elif channel == '+<x>':
            phase_channel = 'BrXPhase'
        elif channel == '-<x>':
            phase_channel = 'BrMinXPhase'
        elif channel == '+<y>':
            phase_channel = 'BrYPhase'
        elif channel == '-<y>':
            phase_channel = 'BrMinYPhase'

        return self.api.hidden[phase_channel].value

    def set_phase(self, channel: str, value: float) -> float:
        if channel == 'Main':
            phase_channel = 'SignalPhase'
        elif channel == '+<x>':
            phase_channel = 'BrXPhase'
        elif channel == '-<x>':
            phase_channel = 'BrMinXPhase'
        elif channel == '+<y>':
            phase_channel = 'BrYPhase'
        elif channel == '-<y>':
            phase_channel = 'BrMinYPhase'

        self.api.hidden[phase_channel].value = value

        return self.get_attenuator(channel)

    def get_ELDOR_freq(self) -> float:
        """Gets the freuency of the ELDOR chnannel.

        Returns:
            float: ELDOR frequency in GHz
        """

        return self.api.cur_exp['ELDORFreqMon'].value

    def set_ELDOR_freq(self, value) -> float:
        """ Sets the freuency of the ELDOR chnannel.

        Args:
            value (float): ELDOR frequency in GHz

        Returns:
            float: ELDOR frequency in GHz
        """

        self.api.cur_exp['ELDORFreqMon'].value = value

        return self.get_ELDOR_freq()

    def get_video_gain(self) -> int:
        return self.api.cur_exp['VideoGain'].value

    def set_video_gain(self, value: int) -> int:
        self.api.cur_exp['VideoGain'].value = value
        return self.get_video_gain()


# =============================================================================
