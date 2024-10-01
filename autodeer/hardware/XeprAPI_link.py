import numpy as np
import time
import os
import sys
import yaml
from autodeer.dataset import create_dataset_from_axes, create_dataset_from_sequence
from scipy.optimize import minimize_scalar
import logging
import re

try:
    xepr_path = os.popen("Xepr --apipath").read()[:-1]
    sys.path.append(xepr_path)
    import XeprAPI

except:
    raise ImportError("Xepr >2.9 needs to be installed on the computer")
hw_log = logging.getLogger('interface.Xepr')

# ============================================================================


class XeprAPILink:
    def __init__(self, config_file: str = None) -> None:
        self.Xepr = None
        self.cur_exp = None
        self.hidden = None
        self._tmp_dir = None
        self.XeprCmds = None
        if config_file is not None:
            with open(config_file, mode='r') as file:
                config = yaml.safe_load(file)
                self.config = config
            self.spec_config = self.config["Spectrometer"]
            self.bridge_config = self.spec_config["Bridge"]
            if self.spec_config["Manufacturer"].lower() \
                    != "bruker":
                msg = "Only Bruker Spectrometers are supported with this this"\
                      "API"
                raise ValueError(msg)
            
            if "Freq Cal" in self.bridge_config.keys():
                self.freq_cal = self.bridge_config["Freq Cal"]
            
            self.AWG = self.spec_config["AWG"]
            self.MPFU = self.spec_config["MPFU"]
            self.Hybrid = self.spec_config["Hybrid"]
        else:
            self.spec_config = {}
            self.AWG = False
            self.MPFU = True
            self.Hybrid = False
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

    def _xepr_retry(self, func, *args, **kwargs):

        for i in range(0, 50):
            try:
                return func(*args, **kwargs)
            except Exception as e:
                time.sleep(1)
                continue
        
        raise Exception("Failed to get Xepr parameter.")

    def get_param(self, param_name: str):
        
        return self._xepr_retry(lambda: self.cur_exp[param_name].value)

    def get_hidden_param(self, param_name: str):
        
        return self._xepr_retry(lambda: self.hidden[param_name].value)
    
    def _set_cur_exp_param(self, param_name, value):
        self.cur_exp[param_name].value = value
    
    def _set_hidden_param(self, param_name, value):
        self.hidden[param_name].value = value
    
    def set_param(self, param_name: str, value):
        self._xepr_retry(self._set_cur_exp_param, param_name, value)
        return self.get_param(param_name)

    def set_hidden_param(self, param_name: str, value):
        self._set_hidden_param(param_name, value)
        # self._xepr_retry(self._set_hidden_param, param_name, value)
        return self.get_hidden_param(param_name)

    def find_Xepr(self):
        open_xepr_instances = XeprAPI.getXeprInstances()

        if len(open_xepr_instances) < 1:
            error_msg = "No Open Xepr API Instances, please open one by" +\
                        ":\n \"processing\" -> \"XeprAPI\" -> \"Enable " +\
                        "Xepr API\""
            hw_log.error(error_msg)
            raise RuntimeError(error_msg)
            
        if self.Xepr is None:
            try:
                self.Xepr_local = XeprAPI.Xepr()
            except OSError:
                error_msg = 'Xepr API: Could not connect to any Xepr instance.'
                hw_log.warning(error_msg)
                print(error_msg)
            except RuntimeError:
                error_msg = "There is already a connection from Xepr to" +\
                            "this python Kernel.\n Please use the correct " +\
                            "python object or restart the kernel "
                hw_log.warning(error_msg)
                raise RuntimeError(error_msg)
            except:
                error_msg = "Can't connect to Xepr: Please check Xepr is " +\
                            "running and open to API"
                hw_log.warning(error_msg)
                raise RuntimeError(error_msg)
            else:
                self._set_Xepr_global(self.Xepr_local)
        else:
            hw_log.debug('Already Connected to Xepr!')
            print('Already Connected to Xepr!')

    def _set_cur_exp_global(self, cur_exp_inst):
        self.cur_exp = cur_exp_inst

    def _get_cur_exp_global(self):
        if self.cur_exp is not None:
            return self.cur_exp
        else:
            error_msg = "Can't find current experiment"
            hw_log.error(error_msg)
            raise RuntimeError(error_msg)

    def find_cur_exp(self):
        """
        Try and finds the current experiment
        """
        if self.Xepr is None:
            error_msg = "Please connect API to Xepr"
            hw_log.error(error_msg)
            raise RuntimeError(error_msg)

        try:
            self.cur_exp = self.Xepr.XeprExperiment()
        except:
            warn_msg = \
                "Can't find the current experiment. Attempting to load it"
            hw_log.warning(warn_msg)
            print(warn_msg)
            self.Xepr.XeprCmds.aqExpSelect("Experiment")
            try:
                self.cur_exp = self.Xepr.XeprExperiment()
            except:
                error_msg = "Can't find the current experiment. Please create"\
                             + "an experiment with the name 'Experiment'"
                hw_log.error(error_msg)
                RuntimeError(error_msg)
            else:
                hw_log.debug("Experiment found")
                print("Experiment found")
        self._set_cur_exp_global(self.cur_exp)
        return self.cur_exp

    def find_hidden(self):
        if self.Xepr is not None:
            self.hidden = self.Xepr.XeprExperiment("AcqHidden")

    def is_exp_running(self):
        return self._xepr_retry(lambda: self.cur_exp.isRunning)

    def acquire_dataset(self, sequence = None):
        """
        This function acquire the dataset, this work both for a running 
        experiment or once it has finished.
        """
        dataclass = self._xepr_retry(lambda: self.Xepr.XeprDataset())
        size = dataclass.size
        data_dim = len(size)
        data = dataclass.O
        params = {
            "nAvgs": self.get_param("recorder.NbScansDone"),
            "reptime": self.get_param("ftEPR.ShotRepTime"),
            "shots": self.get_param("ftEPR.ShotsPLoop"),
            "B": self.get_field(),
            "LO": self.get_counterfreq(),
            }

        if sequence is None:
            default_labels = ['X','Y','Z','T']
            dims = default_labels[:data_dim]
            labels = []

            for i in range(data_dim):
                ax_label = default_labels[i]
                axis_string = self.cur_exp[f'recorder.Abs{i+1}Data'].aqGetParUnits()
                if "'" in axis_string:
                    axis_string = axis_string.replace("'", "")
                if axis_string == 'G':
                    labels.append('B')
                elif axis_string == 'ns':
                    labels.append('t')
                else:
                    labels.append(None)

            if data_dim == 1:
                # We have a 1D dataset
                t = dataclass.X
                hw_log.debug('Acquired Dataset')
                return create_dataset_from_axes(data, t, params,labels)
            elif data_dim == 2:
                # we have a 2D dataset
                t1 = dataclass.X
                t2 = dataclass.Y
                hw_log.debug('Acquired Dataset')
                return create_dataset_from_axes(data, [t1,t2], params,labels)
            
        else:
            dset = create_dataset_from_sequence(data, sequence)
            for par in params:
                dset.attrs[par] = params[par]
            return dset


    def acquire_scan(self,sequence = None):
        """
        This script detects the end of the scan and acquires the data set. 
        This requires that the experiment is still running, or recently 
        finished. Once it has been saved this will no longer work.
        """
        if self.is_exp_running():
            
            self.pause_exp()
            while self.is_exp_running():
                time.sleep(1)
            dataset = self.acquire_dataset(sequence)

            self.rerun_exp()
            time.sleep(0.5)
            return dataset
        else:
            return self.acquire_dataset(sequence)

    def acquire_scan_at(self, scan_num: int,sequence = None):
        """
        This script acquires the scan after a specific number of scans
        """
        x_length = int(self.get_param("XSpecRes"))
        time_per_point = self.get_param("ShotRepTime") * 1e-6 \
            * self.get_param("ShotsPLoop") * 2
        while self.get_param("NbScansDone") != scan_num:
            time.sleep(time_per_point * x_length / 2)
        return self.acquire_scan(sequence)

    def acquire_scan_2d(self,sequence = None):
        """
        This function acquires the dataset after a full 2D scan.
        This is done by identifying the number of scan steps per sweep and 
        acquiring the data on that scan.
        This requires that the experiment has not been saved. 
        """
        if self.is_exp_running():
            total_num_scan = self.get_param("NbScansToDo")
            total_num_sweeps = self.get_param("SweepsPExp")
            scans_per_sweep = total_num_scan/total_num_sweeps

            if not scans_per_sweep.is_integer():
                raise RuntimeError("There is a non integer number of scans per"
                                   " sweep")

            current_scan = self.get_param("NbScansDone")
            current_sweep = np.floor(current_scan/scans_per_sweep)
            next_scan_target = (current_sweep + 1) * scans_per_sweep

            return self.acquire_scan_at(next_scan_target,sequence)
        else:
            return self.acquire_scan(sequence)

    def set_PulseSpel_var(self, variable: str, value: int):
        """
        This can be used to change any pulse spell variable whilst the 
        experiment is running. These changes take effect at the beginning of 
        the next scan
        """

        if self.is_exp_running():
            raise RuntimeWarning('Variable can not be set when experiment is running')
            hw_log.warning(f"Variable {variable} not set as experiment is running.")

        if re.search(r"[pd]\d+", variable) is not None:
            if value % self.bridge_config["Pulse dt"] != 0:
                msg = "Pulse/delay parameter must be of same resolution as" \
                    + "patternjet/AWG. Resolution is" \
                    + f" {self.bridge_config['Pulse dt']} ns"

                raise ValueError(msg)

        hw_log.debug(f'Set pulse spell var {str(variable)} to {int(value)}')
        self.set_param("ftEpr.PlsSPELSetVar", str(variable) + " = " + str(int(value)))

    def get_ReplaceMode(self):
        return self.get_param("ftEpr.ReplaceMode") == "On"

    def set_ReplaceMode(self, state=False):
        init_state = self.get_ReplaceMode()
        if state:
            value = 'On'
            hw_log.warning('Replace mode turned on')
            print('DANGER: Replace Mode turned ON!')
        else:
            value = 'Off'
            if init_state:
                hw_log.info('Replace mode turned off')

        self.set_param("ftEpr.ReplaceMode", value)
        return self.get_ReplaceMode()

    def get_PhaseCycle(self):
        return self.get_param("PCycleOn")

    def set_PhaseCycle(self, state=True):
        init_state = self.get_PhaseCycle()
        if state != init_state:
            hw_log.info(f"On-Board Phase Cycling set to {state}")
            self.set_param("PCycleOn", state)  
        return self.get_PhaseCycle()

    def get_PulseSpel_exp_filename(self):
        return os.path.basename(self.get_param("ftEPR.PlsSPELPrgPaF"))

    def get_PulseSpel_exp_filepath(self):
        self.save_PulseSpel_exp()
        return self.get_param("ftEPR.PlsSPELPrgPaF")

    def set_PulseSpel_exp_filepath(self, fullpath):
        self.save_PulseSpel_exp()
        hw_log.info(f"Set PulseSpel experiment file to {fullpath}")
        self.Xepr.XeprCmds.aqPgLoad(fullpath)

    def get_PulseSpel_def_filename(self):
        return os.path.basename(self.get_param("ftEPR.PlsSPELGlbPaF")) 

    def get_PulseSpel_def_filepath(self):
        return self.get_param("ftEPR.PlsSPELGlbPaF")

    def set_PulseSpel_def_filepath(self, fullpath):
        self.save_PulseSpel_def()
        hw_log.info(f"Set PulseSpel definition file to {fullpath}")
        self.Xepr.XeprCmds.aqPgDefLoad(fullpath)

    def get_PulseSpel_phase_cycling(self):
        return self.get_param("ftEPR.PlsSPELLISTSlct")

    def set_PulseSpel_phase_cycling(self, name):
        hw_log.info(f"Set PulseSpel phase cycle to {name}")
        self.set_param("ftEPR.PlsSPELLISTSlct", name)

    def get_PulseSpel_experiment(self):
        return self.get_param("ftEPR.PlsSPELEXPSlct")

    def set_PulseSpel_experiment(self, name):
        hw_log.info(f"Set PulseSpel experiment to {name}")
        self.set_param("ftEPR.PlsSPELEXPSlct",name)

    def save_PulseSpel_exp(self, name: str = None):
        if name is None:
            # Save as a temp file
            if self._tmp_dir is None:
                try:
                    os.mkdir('/tmp/autoDeer/')
                    self._tmp_dir = '/tmp/autoDeer/'
                except FileExistsError:
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
                except FileExistsError:
                    self._tmp_dir = '/tmp/autoDeer/'
                    print("temp directory already exists")
            timestr = time.strftime("%Y%m%d-%H%M%S")
            path = self._tmp_dir + "pulsespel_" + timestr + ".def"
        else:
            path = name
        hw_log.debug(f'Saved Pulse Spel definition to: {path}')
        self.Xepr.XeprCmds.aqPgDefSaveAs(path)

    def get_Acquisition_mode(self):
        return self.get_param("ftEPR.FTAcqModeSlct")

    def set_Acquisition_mode(self, mode: int):
        """mode=0: Run from tables, mode=1: Run from Pulse Spel,
        mode=2:Read transient, mode=3:Start Transient"""
        if mode == 0:
            self.set_param("ftEPR.FTAcqModeSlct",'Run from Tables')
        elif mode == 1:
            self.set_param("ftEPR.FTAcqModeSlct",'Run from PulseSPEL')
        elif mode == 2:
            self.set_param("ftEPR.FTAcqModeSlct",'Read Transient')
        elif mode == 3:
            self.set_param("ftEPR.FTAcqModeSlct",'Start Transient')
        else:
            print('Acquisition Mode not changed. Input error.')
            return 0
        hw_log.debug(f"Acquistiton mode set to: {mode}")

    def compile_PulseSpel_prg(self):
        self.Xepr.XeprCmds.aqPgShowPrg()
        self.Xepr.XeprCmds.aqPgCompile()
        hw_log.debug("PulseSpel Exp file compiled")
        time.sleep(0.5)
        pass

    def compile_PulseSpel_def(self):
        self.Xepr.XeprCmds.aqPgShowDef()
        self.Xepr.XeprCmds.aqPgCompile()
        hw_log.debug("PulseSpel Exp file compiled")
        time.sleep(0.5)
        pass

    def run_exp(self):
        """
        Runs the current experiment.
        """
        self.cur_exp.aqExpRun()
        hw_log.debug('Experiment started')
        time.sleep(5)  # Required
        pass

    def pause_exp(self):
        """
        Pauses the current experiment.
        """
        self.cur_exp.aqExpPause()
        hw_log.debug('Experiment Paused')
        pass

    def rerun_exp(self):
        """
        Re-runs the current experiment. Only for use after a pause commands
        """
        self.cur_exp.aqExpRun()
        hw_log.debug('Experiment rerun')
        pass

    def stop_exp(self):
        """
        Stops the current experiment.
        """
        self.cur_exp.aqExpStop()
        hw_log.debug('Experiment stopped')
        pass

    def abort_exp(self):
        """
        Aborts the current experiment.
        """
        self.cur_exp.aqExpAbort()
        hw_log.debug('Experiment aborted')
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
        xepr_file_limit = 120
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

    def get_field(self) -> int:
        """Gets the magnetic field.

        Returns
        -------
        int
            Field position in Gauss
        """
        return self.get_param('CenterField')

    def set_field(self, val: int, hold: bool = True) -> int:
        """Sets the magentic field.

        Parameters
        ----------
        val : int
            Field position in Gauss
        hold : bool, optional
            Wait until Xepr belives the field is stable, by default True

        Returns
        -------
        int
            Field position in Gauss
        """
        self.set_param('CenterField', np.around(val,3))
        time.sleep(2)  # Always wait 2s after a field change
        hw_log.info(f'Field position set to {val} G')
        if hold is True:
            while self.get_param('FieldWait') is True:
                time.sleep(0.5)
        return self.get_field()

    def get_counterfreq(self) -> float:
        """Gets the frequency counter value.

        Note: There is no equivalent set command 

        Returns
        -------
        float
            The counter frequency in GHz
        """
        return self.get_param('FrequencyMon')

    def set_sweep_width(self, val: int) -> int:
        """Sets the field sweep width

        Parameters
        ----------
        val : int
            Field sweep width in Gauss

        Returns
        -------
        int
            Field sweep width in Gauss
        """
        self.set_param('SweepWidth',val)
        hw_log.debug('Field sweep width set to {val:.6g} G')
        return self.get_sweep_width()
    
    def get_sweep_width(self) -> int:
        """Gets the field sweep width

        Returns
        -------
        int
            Field sweep width in Gauss.
        """
        return self.get_param('SweepWidth')

    def set_freq(self, val: float, pol: list = None,
                 precision: bool = False) -> float:
        """Sets the current bridge frequency. The stepper motor value is 
        calculated through a conversion. 

        Parameters
        ----------
        val : float
            The frequency in GHz
        pol : list, optional
            Conversion polynomal coefficents for freq->Step, by default None

        Returns
        -------
        float
            Frequency in GHz
        """
        if "Min Freq" in self.bridge_config.keys():
            if val < self.bridge_config["Min Freq"]:
                raise RuntimeError("Set Frequency is too low!")

        if "Max Freq" in self.bridge_config.keys():
            if val > self.bridge_config["Max Freq"]:
                raise RuntimeError("Set Frequency is too high!")
                
        if "Digital Source" in self.bridge_config.keys() and self.bridge_config["Digital Source"]:
            if self.bridge_config["Digital Source"]:
                self.set_hidden_param('FineFreq',val)
                hw_log.info(f'Bridge Frequency set to {val}Fset')
        else:
            if pol is None:
                if hasattr(self, "freq_cal"):
                    pol = self.freq_cal
                else:
                    pol = [-67652.70, 2050.203]

            pol_func = np.polynomial.polynomial.Polynomial(pol)
            pos = round(pol_func(val))
            if precision:
                self.set_hidden_param('Frequency',pos)
                bounds = [pos-50, pos+50]

                def objective(x):
                    self.set_hidden_param('Frequency',x)
                    return self.get_counterfreq()

                output = minimize_scalar(
                    objective, method='bounded', bounds=bounds,
                    options={'xatol': 1e-3, 'maxiter': 30})
                pos = round(output.x)
                self.set_hidden_param('Frequency',pos)
            else:
                self.set_hidden_param('Frequency',pos)
    
            hw_log.info(f'Bridge Frequency set to {val} at position {pos}')

        return self.get_hidden_param('Frequency')

    def get_freq(self, pol: list = None) -> float:
        """Gets the current bridge frequency. The frequency is calculated
        through a conversion. 

        Parameters
        ----------
        pol : list, optional
            Conversion polynomal coefficents for freq->Step, by default None

        Returns
        -------
        float
            Frequency in GHz
        """
        if "Digital Source" in self.bridge_config.keys():
            if self.bridge_config["Digital Source"]:
                return self.get_hidden_param('FineFreq')
        else:
            if pol is None:
                if hasattr(self, "freq_cal"):
                    pol = self.freq_cal
                else:
                    pol = [-67652.70, 2050.203]
            inv_func = np.polynomial.polynomial.Polynomial(
                [-pol[0]/pol[1], 1/(pol[1])])

            return inv_func(self.hidden['Frequency'].value)

    def get_spec_config(self) -> str:
        """get_spec_config Gets the name of the current spectrometer
        configuration file.

        Returns
        -------
        str
            Returns the name of the current spectrometer configuration file
        """

        return self.get_hidden_param('PlsPrgCalDbName')

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

        if atten_channel == 'ELDORAtt':
            return self.get_param(atten_channel)
        else:
            return self.get_hidden_param(atten_channel)

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
        
        if atten_channel == 'ELDORAtt':
            self.set_param(atten_channel, value)
        else:
            self.set_hidden_param(atten_channel,value)
        hw_log.debug(f"{channel} Attenuator set to {value} ")

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

        return self.get_hidden_param(phase_channel)

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

        self.set_hidden_param(phase_channel,value)
        hw_log.debug(f"{channel} phase control set to {value} ")

        return self.get_attenuator(channel)

    def get_ELDOR_freq(self) -> float:
        """Gets the freuency of the ELDOR chnannel.

        Returns:
            float: ELDOR frequency in GHz
        """

        return self.get_param('FrequencyA')

    def set_ELDOR_freq(self, value: float) -> float:
        """ Sets the freuency of the ELDOR chnannel.

        Args:
            value (float): ELDOR frequency in GHz

        Returns:
            float: ELDOR frequency in GHz
        """

        if "Min Freq" in self.bridge_config.keys():
            if value < self.bridge_config["Min Freq"]:
                raise RuntimeError("Set ELDOR Frequency is too low!")

        if "Max Freq" in self.bridge_config.keys():
            if value > self.bridge_config["Max Freq"]:
                raise RuntimeError("Set ELDOR Frequency is too high!")

        self.set_param('FrequencyA', value)
        hw_log.debug(f"ELDOR frequency set to {value}")
        return self.get_ELDOR_freq()

    def get_video_gain(self) -> int:
        """Get the video gain in dB

        Returns
        -------
        int
            Video gain in dB
        """
        return self.get_param('VideoGain')
    
    def get_video_gain_step(self) -> int:

        return self.cur_exp['VideoGain'].aqGetParFineSteps()

    def set_video_gain(self, value: int) -> int:
        """Set the video gain in dB

        Parameters
        ----------
        value : int
            Video gain to be set in dB.

        Returns
        -------
        int
            The adjusted video gain
        """
        self.set_param('VideoGain', value)
        hw_log.debug(f"Video gain set to {value}")
        return self.get_video_gain()

    def get_config_file(self) -> str:
        """Get the name of the spectrometer configuration file

        Returns
        -------
        str
            The name of the spectrometer configuration file
        """
        return self.get_hidden_param('PlsPrgCalDbName')

    def set_config_file(self, config_file: str) -> str:
        """Sets the name of the spectrometer configuration file. If the current
        configuration file is already correct, nothing changes.
       
        Parameters
        ----------
        config_file : str
            The new configuration file name.

        Returns
        -------
        str
            The name of the spectrometer configuration file
        """
        if self.get_config_file() != config_file:
            self.set_hidden_param('PlsPrgCalDbName',config_file)
            self.hidden['PlsPrgCalDbLoad']
            self.hidden['ApplyCfg']
            hw_log.debug("Changed config file")
        else:
            hw_log.debug("Not changed config file")

        return self.get_config_file()

    def get_MW_amp(self) -> bool:
        """Gets the current setting of the microwave amplifier

        Returns
        -------
        bool
            current setting of the microwave amplifier
        """
        if self.hidden['MWGain'].value == 'On':
            return 1
        elif self.hidden['MWGain'].value == 'Off':
            return 0
        else:
            return self.hidden['MWGain']

    def set_MW_amp(self, value: bool) -> bool:
        """Sets the current setting of the microwave amplifier

        Parameters
        ----------
        value : bool
            New setting of the microwave amplifier

        Returns
        -------
        bool
            Setting of the microwave amplifier
        """
        if value:
            self.set_hidden_param('MWGain','On')
        else:
            self.set_hidden_param('MWGain','Off')
        hw_log.debug("Microwave amplifier set to {value}")
        return self.get_MW_amp()

    def get_video_bandwidth(self) -> float:
        """Gets the detection video bandwidth.

        Returns
        -------
        float
            Video bandwidth in MHz
        """
        return self.get_video_bandwidth()

    def set_video_bandwidth(self, value: float) -> float:
        """Sets the detection video bandwidth.

        Parameters
        ----------
        value : float
            Video bandwidth in MHz, options = [200,20].

        Returns
        -------
        float
            Video bandwidth in MHz
        """
        options = [200, 20]

        if round(value) not in options:
            raise ValueError("Video Bandwidth can only be 200MHz or 20MHz")

        self.hidden['VideoBW'] = value

        return self.get_video_bandwidth()
# =============================================================================
