from distutils.command.config import config
from multiprocessing.connection import wait
import numpy as np
import time
import os,sys
from numpy.core.fromnumeric import cumprod;
import XeprAPI


import logging

hw_log = logging.getLogger('hardware.Xepr')

# def connect_Xepr():
#     try:

#         Xepr = XeprAPI.Xepr()
#         return Xepr

hardware_meta = {# This dictionary should be moved into a config file eventually
    "Type":             "Complete Spectrometer",
    "Manufactuer":      "Bruker",
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
        self.meta = hardware_meta # This will become more neuanced eventually. Currentlty this is sufficent.
        pass

    def set_Xepr_global(self,Xepr_inst): 
        self.Xepr = Xepr_inst
        self.XeprCmds = self.Xepr.XeprCmds

    def get_Xepr_global(self):
        if self.Xepr != None:
            return self.Xepr
        else:
            raise RuntimeError("Can't find XEPR instance")

    def find_Xepr(self):
        open_xepr_instances = XeprAPI.getXeprInstances()

        if len(open_xepr_instances) < 1:
            raise RuntimeError("No Open Xepr API Instances, please open one by:\n"+\
                "\"processing\" -> \"XeprAPI\" -> \"Enable Xepr API\"")

        if self.Xepr == None:
            try:
                self.Xepr_local = XeprAPI.Xepr()
            except OSError:
                print('Xepr API: Could not connect to any Xepr instance.')
            except RuntimeError:
                raise RuntimeError("There is already a connection from Xepr to this python Kernal.\n" \
                    "Please use the correct python object or restart the kernal ")
            except:
                RuntimeError("Can't connect to Xepr: Please check Xepr is running and open to API")
            else:
                self.set_Xepr_global(self.Xepr_local)
        else:
            print('Already Connected to Xepr!')
            
    def set_cur_exp_global(self,cur_exp_inst):
        self.cur_exp = cur_exp_inst

    def get_cur_exp_global(self):
        if self.cur_exp != None:
            return self.cur_exp
        else:
            raise RuntimeError("Can't find current experiment")
            
    def find_cur_exp(self):
        """
        Trys and finds the current experiment
        """
        if self.Xepr == None:
            raise RuntimeError("Please connect API to Xepr")
        
        try:
            self.cur_exp = self.Xepr.XeprExperiment()
        except:
            print("Can't find the current experiment. Attempting to load it")
            self.Xepr.XeprCmds.aqExpSelect("Experiment")
            try:
                self.cur_exp = self.Xepr.XeprExperiment()
            except:
                RuntimeError("Can't find the current experiment. Please create an experiment with name 'Experiment'")
            else:
                print("Experiment found")
        self.set_cur_exp_global(self.cur_exp)
        return self.cur_exp

    def find_hidden(self):
        if self.Xepr != None:
            self.hidden = self.Xepr.XeprExperiment("AcqHidden")

    def is_exp_running(self):
        return self.cur_exp.isRunning

    def acquire_dataset(self): 
        """
        This function acquire the dataset, this work both for a running experiment or once it has finished.
        """
        dataclass = self.Xepr.XeprDataset() # This function returns the currently view dataset
        size = dataclass.size # This needs to be checked to see what form this precisely comes as
        data_dim = len(size)
        data = dataclass.O
        if data_dim == 1:
            # We have a 1D dataset
            t = dataclass.X
            hw_log.debug('Acquired Dataset')
            return dataset(t,data,self.cur_exp)
        elif data_dim == 2:
            # we have a 2D dataset
            t1 = dataclass.X
            t2 = dataclass.Y
            hw_log.debug('Acquired Dataset')
            return dataset([t1,t2],data,self.cur_exp)

    def acquire_scan(self):
        """
        This script detects the end of the scan and acquires the data set. This requries that the experiment is still running, or recently finished. Once it has been saved this will no longer work.
        """
        if self.is_exp_running():
        
            current_scan = self.cur_exp.getParam("NbScansDone").value
            x_length = int(self.cur_exp.getParam("XSpecRes").value)
            time_per_point = self.cur_exp.getParam("ShotRepTime").value *1e-6*self.cur_exp.getParam("ShotsPLoop").value*2
            trace = np.zeros(x_length, dtype = np.complex64)
            while self.cur_exp.getParam("NbScansDone").value == current_scan:
                time.sleep(time_per_point)
            time.sleep(time_per_point)
            return self.acquire_dataset()
        else:
            return self.acquire_dataset()

    def acquire_scan_at(self,scan_num:int):
        """
        This script acquires the scan after a spacific number of scans
        """
        x_length = int(self.cur_exp.getParam("XSpecRes").value)
        time_per_point = self.cur_exp.getParam("ShotRepTime").value *1e-6*self.cur_exp.getParam("ShotsPLoop").value*2
        while self.cur_exp.getParam("NbScansDone").value != scan_num:
            time.sleep(time_per_point*x_length/2)
        return self.acquire_scan()

    def acquire_scan_2d(self):
        """
        This function acquries the dataset after a full 2D scan.
        This is done by identfying the number of scansteps per sweep and acquring the data on that scan.
        This requires that the experiment has not been saved. 
        """
        if self.is_exp_running():
            total_num_scan = self.cur_exp.getParam("NbScansToDo").value
            total_num_sweeps = self.cur_exp.getParam("SweepsPExp").value
            scans_per_sweep = total_num_scan/total_num_sweeps

            if not scans_per_sweep.is_integer():
                raise RuntimeError('There is a non integer number of scans per sweep')

            current_scan = self.cur_exp.getParam("NbScansDone").value
            current_sweep = np.floor(current_scan/scans_per_sweep)
            next_scan_target = (current_sweep + 1) * scans_per_sweep

            return self.acquire_scan_at(next_scan_target)
        else:
            return self.acquire_scan()


    def set_PulseSpel_var(self,variable:str,value:int):
        """
        This can be used to change any pulse spel variable whilst the experiment is running. These changes take effect at the begining of the next scan
        """
        hw_log.debug(f'Set pulse spel var {str(variable)} to {int(value)}')
        self.cur_exp["ftEpr.PlsSPELSetVar"].value = str(variable) + " = "+ str(int(value))
    
    def set_ReplaceMode(self,state=False):
        if state:
            value = 'On'
            hw_log.warning('Replace mode turned on')
            print('DANGER: Replace Mode turned ON!')
        else:
            value = 'Off'
            hw_log.info('Replace mode turned off')

        self.cur_exp["ftEpr.ReplaceMode"].value = value

    def set_PhaseCycle(self,state=True):
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
        return self.cur_exp["ftEPR.PlsSPELPrgPaF"].value

    def set_PulseSpel_exp_filepath(self,fullpath):
        self.save_PulseSpel_exp()
        self.Xepr.XeprCmds.aqPgLoad(fullpath)

    def get_PulseSpel_def_filename(self):
        return os.path.basename(self.cur_exp["ftEPR.PlsSPELGlbPaF"].value) 

    def get_PulseSpel_def_filenpath(self):
        return self.cur_exp["ftEPR.PlsSPELGlbPaF"].value

    def set_PulseSpel_def_filepath(self,fullpath):
        self.save_PulseSpel_def()
        self.Xepr.XeprCmds.aqPgDefLoad(fullpath)

    def get_PulseSpel_phase_cycling(self):
        return self.cur_exp["ftEPR.PlsSPELLISTSlct"].value

    def set_PulseSpel_phase_cycling(self,name):
        self.cur_exp["ftEPR.PlsSPELLISTSlct"].value = name
        
        if self.cur_exp["ftEPR.PlsSPELLISTSlct"].value != name:
            print("WARNING: Phase cycling did not change")
            return 0
        else:
            return 1

    def get_PulseSpel_experiment(self):
        return self.cur_exp["ftEPR.PlsSPELEXPSlct"].value
    
    def set_PulseSpel_experiment(self,name):
        self.cur_exp["ftEPR.PlsSPELEXPSlct"].value = name
        
        if self.cur_exp["ftEPR.PlsSPELEXPSlct"].value != name:
            print("WARNING: Pulse Spel Experiment did not change")
            return 0
        else:
            return 1
    

    def save_PulseSpel_exp(self,name=None):
        if name==None:
            # Save as a temp file
            if self._tmp_dir == None:
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

    def save_PulseSpel_def(self,name=None):
        if name == None:
            # Save as a temp file
            if self._tmp_dir == None:
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
        hw_log.debug(f'Saved Pulse Spel defintion to: {path}')
        self.Xepr.XeprCmds.aqPgDefSaveAs(path)



    def get_Acquistion_mode(self):
        return self.cur_exp["ftEPR.FTAcqModeSlct"].value    
    
    def set_Acquistion_mode(self, mode:int):
        """mode=0: Run from tabels, mode=1: Run from Pulse Spel, mode=2:Read transient, mode=3:Start Transient"""
        if mode == 0:
            self.cur_exp["ftEPR.FTAcqModeSlct"].value = 'Run from Tabels'
        elif mode == 1:
            self.cur_exp["ftEPR.FTAcqModeSlct"].value = 'Run from PulseSPEL'
        elif mode == 2:
            self.cur_exp["ftEPR.FTAcqModeSlct"].value = 'Read Transient'
        elif mode == 3:
            self.cur_exp["ftEPR.FTAcqModeSlct"].value = 'Start Transient'
        else:
            print('Acqusiton Mode not changed. Input error.')
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
        time.sleep(5)
        pass

    def stop_exp(self):
        self.cur_exp.aqExpStop()
        hw_log.info('Experiment stopped')
        pass

    def abort_exp(self):
        self.cur_exp.aqExpAbort()
        hw_log.info('Experiment aborted')
        pass

    def xepr_save(self,path,title=None):
        # Saves the current viewpoint to either the specified file in the working directory or to the filepath.
        # This is a bruker save function
        # Taken from Custom Xepr
        directory, basename = os.path.split(path)
        if not title:
            title = basename
        if not os.path.exists(directory):
            os.makedirs(directory)
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

    def set_field(self,val:int,hold:bool=True) -> int:
        """ This sets the central field"""
        self.cur_exp['CenterField'].value = val
        time.sleep(2) #Always wait 2s after a field change
        hw_log.info(f'Field position set to {val} G')
        if hold == True:
            while self.cur_exp['FieldWait'] == True:
                time.sleep(0.5)
        return self.cur_exp['CenterField'].value

    def get_counterfreq(self) -> float:
        """ This returns the current freq counter"""
        return self.cur_exp['FrequencyMon'].value

    def set_sweep_width(self,val:int) -> int:
        self.cur_exp['SweepWidth'].value = val
        hw_log.info('Field sweep width set to {val} G')
        return self.cur_exp['SweepWidth'].value 
    
    def get_sweep_width(self) -> int:
        return self.cur_exp['SweepWidth'].value

    def set_freq(self,val:np.float128) -> float:
        """ This sets bridge frequency, and works through a polynomial approximation. This might need to be adjusted
        for different spectrometers"""
        f_pol = [1.420009750632201e4,
            -5.118516287710228e3,
            2.092103562165744e2,
            -2.034307248428457,
            0,0]

        pol_func = np.polynomial.polynomial.Polynomial(f_pol)
        pos = round(pol_func(val))
        self.hidden['Frequency'].value = pos
        hw_log.info(f'Frequency set to {val} at position {pos}')
        return self.hidden['Frequency'].value

    def get_freq(self) -> float:
        """ This returns the current bridge frequency"""
        return self.hidden['Frequency'].value

    def get_spec_config(self) -> str:
        """get_spec_config Gets the name of the current spectrometer configuration file.

        Returns
        -------
        str
            Returns the name of the current spectrometer configuration file
        """

        return self.hidden['PlsPrgCalDbName'].value

    def set_spec_config(self,name:str='Normal') -> str:
        """set_spec_config Sets the name of the current spectrometer configuration file

        Parameters
        ----------
        name : str, optional
            The file name of config file. Normal and AWG, shortcut to the standard types, by default 'Normal'

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
        



    


## Section on phase control

class phase:
    """
    A class for the control of phases.
    Options: "cwBridge.SignalPhase", "ftBridge.BrXPhase", "ftBridge.BrYPhase", "ftBridge.BrMinXPhase", "ftBridge.BrMinYPhase"
    """
    def __init__(self,name:str,hidden):
        self.name = name
        self.hidden = hidden
        self.max = self.hidden[self.name].aqGetParMaxValue
        self.min = self.hidden[self.name].aqGetParMinValue
        self.course_step = self.hidden[self.name].aqGetParCoarseSteps
        self.fine_step = self.hidden[self.name].aqGetParFineSteps
    
    def get_value(self):
        return self.hidden[self.name].value
    
    def set_value(self,new_value):
        self.hidden[self.name.value] = new_value
        return self.hidden[self.name].value

class attenuator:
    """
    A class for the control of both stepped and variable attenuators
    Name(str) - The Xepr code for the attenuator. Options: "ftBridge.BrXAmp", "ftBridge.BrYAmp", "ftBridge.BrMinXAmp", "ftBridge.BrMinYAmp"
    """
    def __init__(self,name:str,hidden):
        self.name = name
        self.hidden = hidden
        self.max = self.hidden[self.name].aqGetParMaxValue
        self.min = self.hidden[self.name].aqGetParMinValue
        self.course_step = self.hidden[self.name].aqGetParCoarseSteps
        self.fine_step = self.hidden[self.name].aqGetParFineSteps
    
    def get_value(self):
        return self.hidden[self.name].value
    
    def set_value(self,new_value):
        self.hidden[self.name.value] = new_value
        return self.hidden[self.name].value         
               
class dataset:
    """ 
    This is the dataset object for openEPR, it not only contains the raw data but also some 
    dataset specific metadata such as current number of scans
    """

    def __init__(self,time,data,cur_exp=None) -> None:
        self.time = time
        self.data = data
        self.dim:int = int(len(time))
        self.size = np.shape(data)
        if cur_exp != None:
            self.scans_done = cur_exp.getParam("NbScansDone").value
            self.scans_todo = cur_exp.getParam("NbScansToDo").value
            self.shrt = cur_exp.getParam("ShotRepTime").value
            self.shot_p_point = cur_exp.getParam("ShotsPLoop").value
        pass


