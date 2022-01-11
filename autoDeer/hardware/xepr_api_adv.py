from posixpath import expanduser
import numpy as np
import time
import os,sys
from numpy.core.fromnumeric import cumprod;
import XeprAPI


# def connect_Xepr():
#     try:
#         Xepr = XeprAPI.Xepr()
#         return Xepr
class xepr_api:
    def __init__(self) -> None:
        self.Xepr = None
        self.cur_exp = None
        self.hidden = None
        pass

    def set_Xepr_global(self,Xepr_inst): 
        self.Xepr = Xepr_inst

    def get_Xepr_global(self):
        if self.Xepr != None:
            return self.Xepr
        else:
            raise RuntimeError("Can't find XEPR instance")

    def find_Xepr(self):
        try:
            self.Xepr_local = XeprAPI.Xepr()
        except OSError:
            print('Xepr already connected')
        except:
            RuntimeError("Can't connect to Xepr: Please check Xepr is running and open to API")
        else:
            self.set_Xepr_global(self.Xepr_local)
            
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
        dataset = self.Xepr.XeprDataset() # This function returns the currently view dataset
        size = dataset.size # This needs to be checked to see what form this precisely comes as
        dataset_dim = len(size)
        data = dataset.O
        if dataset_dim == 1:
            # We have a 1D dataset
            t = dataset.X
            return t,data
        elif dataset_dim == 2:
            # we have a 2D dataset
            t1 = dataset.X
            t2 = dataset.Y
            return t1,t2,data

    def acquire_scan(self):
        """
        This script detects the end of the scan and acquires the data set. This requries that the experiment is still running, or recently finished. Once it has been saved this will no longer work.
        """
        current_scan = self.cur_exp.getParam("NbScansDone").value
        x_length = int(self.cur_exp.getParam("XSpecRes").value)
        time_per_point = self.cur_exp.getParam("ShotRepTime").value *1e-6*self.cur_exp.getParam("ShotsPLoop").value*2
        trace = np.zeros(x_length, dtype = np.complex64)
        while self.cur_exp.getParam("NbScansDone").value == current_scan:
            time.sleep(time_per_point)
        time.sleep(time_per_point)
        # trace = Xepr.().O
        # time_axis = Xepr.XeprDataset().X
        # return time_axis,trace
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
        total_num_scan = self.cur_exp.getParam("NbScansToDo").value
        total_num_sweeps = self.cur_exp.getParam("SweepsPExp").value
        scans_per_sweep = total_num_scan/total_num_sweeps

        if not scans_per_sweep.is_integer():
            raise RuntimeError('There is a non integer number of scans per sweep')
        
        current_scan = self.cur_exp.getParam("NbScansDone").value
        current_sweep = np.floor(current_scan/scans_per_sweep)
        next_scan_target = (current_sweep + 1) * scans_per_sweep

        return self.acquire_scan_at(next_scan_target)


    def set_PulseSpel_var(self,cur_exp,variable:str,value:int):
        """
        This can be used to change any pulse spel variable whilst the experiment is running. These changes take effect at the begining of the next scan
        """
        cur_exp["ftEpr.PlsSPELSetVar"].value = str(variable) + " = "+ str(int(value))
    
    def set_ReplaceMode(self,cur_exp,state=False):
        if state:
            value = 'On'
            print('DANGER: Replace Mode turned ON!')
        else:
            value = 'Off'
        cur_exp["ftEpr.ReplaceMode"].value = value
    
    def get_PulseSpel_exp_filename(self):
        return os.path.basename(self.cur_exp["ftEPR.PlsSPELPrgPaF"].value)

    def get_PulseSpel_exp_filepath(self):
        return self.cur_exp["ftEPR.PlsSPELPrgPaF"].value

    def set_PulseSpel_exp_filepath(self,fullpath):
        self.Xepr.XeprCmds.aqPgLoad(fullpath)

    def get_PulseSpel_def_filename(self):
        return os.path.basename(self.cur_exp["ftEPR.PlsSPELGlbPaF"].value) 

    def get_PulseSpel_def_filenpath(self):
        return self.cur_exp["ftEPR.PlsSPELGlbPaF"].value

    def set_PulseSpel_def_filepath(self,fullpath):
        self.Xepr.XeprCmds.aqPgDefLoad(fullpath)

    def get_PulseSpel_phase_cycling(self):
        return self.cur_exp["ftEPR.PlsSPELLISTSlct"].value

    def set_PulseSpel_phase_cycling(self,cur_exp,name):
        self.cur_exp["ftEPR.PlsSPELLISTSlct"].value = name
        
        if self.cur_exp["ftEPR.PlsSPELLISTSlct"].value != name:
            print("WARNING: Phase cycling did not change")
            return 0
        else:
            return 1

    def get_PulseSpel_experiment(self,cur_exp):
        return cur_exp["ftEPR.PlsSPELEXPSlct"].value
    
    def set_PulseSpel_experiment(self,cur_exp,name):
        cur_exp["ftEPR.PlsSPELEXPSlct"].value = name
        
        if cur_exp["ftEPR.PlsSPELEXPSlct"].value != name:
            print("WARNING: Pulse Spel Experiment did not change")
            return 0
        else:
            return 1

    def get_Acquistion_mode(self,cur_exp):
        return cur_exp["ftEPR.FTAcqModeSlct"].value    
    
    def set_Acquistion_mode(self,cur_exp, mode:int):
        """mode=0: Run from tabels, mode=1: Run from Pulse Spel, mode=2:Read transient, mode=3:Start Transient"""
        if mode == 0:
            cur_exp["ftEPR.FTAcqModeSlct"].value = 'Run from Tabels'
        elif mode == 1:
            cur_exp["ftEPR.FTAcqModeSlct"].value = 'Run from PulseSPEL'
        elif mode == 2:
            cur_exp["ftEPR.FTAcqModeSlct"].value = 'Read Transient'
        elif mode == 3:
            cur_exp["ftEPR.FTAcqModeSlct"].value = 'Start Transient'
        else:
            print('Acqusiton Mode not changed. Input error.')
            return 0
              
    def compile_PulseSpel_prg(self):
        self.Xepr.XeprCmds.aqPgShowPrg()
        self.Xepr.XeprCmds.aqPgCompile()
        time.sleep(0.5)

    def compile_PulseSpel_def(self):
        self.Xepr.XeprCmds.aqPgShowDef()
        self.Xepr.XeprCmds.aqPgCompile()
        time.sleep(0.5)


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
               