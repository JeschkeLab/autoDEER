from posixpath import expanduser
import numpy as np
import time
import os,sys
from numpy.core.fromnumeric import cumprod;
import XeprAPI
import deerlab as dl


# def connect_Xepr():
#     try:
#         Xepr = XeprAPI.Xepr()
#         return Xepr
Xepr = None
cur_exp = None
hidden = None

def set_Xepr_global(Xepr_inst): 
    global Xepr
    Xepr = Xepr_inst

def get_Xepr_global():
    if Xepr != None:
        return Xepr
    else:
        raise RuntimeError("Can't find XEPR instance")

def find_Xepr():
    try:
        Xepr_local = XeprAPI.Xepr()
    except OSError:
        print('Xepr already connected')
    except:
        RuntimeError("Can't connect to Xepr: Please check Xepr is running and open to API")
    else:
        set_Xepr_global(Xepr_local)
        
def set_cur_exp_global(cur_exp_inst):
    global cur_exp
    cur_exp = cur_exp_inst

def get_cur_exp_global():
    if cur_exp != None:
        return cur_exp
    else:
        raise RuntimeError("Can't find current experiment")
        
def find_cur_exp():
    """
    Trys and finds the current experiment
    """
    if Xepr == None:
        raise RuntimeError("Please connect API to Xepr")
    
    try:
        cur_exp = Xepr.XeprExperiment()
    except:
        print("Can't find the current experiment. Attempting to load it")
        Xepr.XeprCmds.aqExpSelect("Experiment")
        try:
            cur_exp = Xepr.XeprExperiment()
        except:
            RuntimeError("Can't find the current experiment. Please create an experiment with name 'Experiment'")
        else:
            print("Experiment found")
    set_cur_exp_global(cur_exp)
    return cur_exp

def find_hidden():
    global hidden
    if Xepr != None:
        hidden = Xepr.XeprExperiment("AcqHidden")

def is_exp_running(): # Untested
    return cur_exp.isRunning

def acquire_dataset(): 
    """
    This function acquire the dataset, this work both for a running experiment or once it has finished.
    """
    dataset = Xepr.XeprDataset() # This function returns the currently view dataset
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

    
def acquire_scan():
    """
    This script detects the end of the scan and acquires the data set. This requries that the experiment is still running, or recently finished. Once it has been saved this will no longer work.
    """
    current_scan = cur_exp.getParam("NbScansDone").value
    x_length = int(cur_exp.getParam("XSpecRes").value)
    time_per_point = cur_exp.getParam("ShotRepTime").value *1e-6*cur_exp.getParam("ShotsPLoop").value*2
    trace = np.zeros(x_length, dtype = np.complex64)
    while cur_exp.getParam("NbScansDone").value == current_scan:
        time.sleep(time_per_point)
    time.sleep(time_per_point)
    # trace = Xepr.().O
    # time_axis = Xepr.XeprDataset().X
    # return time_axis,trace
    return acquire_dataset()
def acquire_scan_at(scan_num:int):
    """
    This script acquires the scan after a spacific number of scans
    """
    x_length = int(cur_exp.getParam("XSpecRes").value)
    time_per_point = cur_exp.getParam("ShotRepTime").value *1e-6*cur_exp.getParam("ShotsPLoop").value*2
    while cur_exp.getParam("NbScansDone").value != scan_num:
        time.sleep(time_per_point*x_length/2)
    return acquire_scan()

def acquire_scan_2d():
    """
    This function acquries the dataset after a full 2D scan.
    This is done by identfying the number of scansteps per sweep and acquring the data on that scan.
    This requires that the experiment has not been saved. 
    """
    total_num_scan = cur_exp.getParam("NbScansToDo").value
    total_num_sweeps = cur_exp.getParam("SweepsPExp").value
    scans_per_sweep = total_num_scan/total_num_sweeps

    if not scans_per_sweep.is_integer():
        raise RuntimeError('There is a non integer number of scans per sweep')
    
    current_scan = cur_exp.getParam("NbScansDone").value
    current_sweep = np.floor(current_scan/scans_per_sweep)
    next_scan_target = (current_sweep + 1) * scans_per_sweep

    return acquire_scan_at(next_scan_target)


def set_PulseSpel_var(cur_exp,variable:str,value:int):
    """
    This can be used to change any pulse spel variable whilst the experiment is running. These changes take effect at the begining of the next scan
    """
    cur_exp["ftEpr.PlsSPELSetVar"].value = str(variable) + " = "+ str(int(value))
    
def set_ReplaceMode(cur_exp,state=False):
    if state:
        value = 'On'
        print('DANGER: Replace Mode turned ON!')
    else:
        value = 'Off'
    cur_exp["ftEpr.ReplaceMode"].value = value
    
def get_PulseSpel_exp_filename():
    return os.path.basename(cur_exp["ftEPR.PlsSPELPrgPaF"].value)

def get_PulseSpel_exp_filepath():
    return cur_exp["ftEPR.PlsSPELPrgPaF"].value

def set_PulseSpel_exp_filepath(fullpath):
    Xepr.XeprCmds.aqPgLoad(fullpath)

def get_PulseSpel_def_filename():
    return os.path.basename(cur_exp["ftEPR.PlsSPELGlbPaF"].value) 

def get_PulseSpel_def_filenpath():
    return cur_exp["ftEPR.PlsSPELGlbPaF"].value

def set_PulseSpel_def_filepath(fullpath):
    Xepr.XeprCmds.aqPgDefLoad(fullpath)

def get_PulseSpel_phase_cycling():
    return cur_exp["ftEPR.PlsSPELLISTSlct"].value

def set_PulseSpel_phase_cycling(cur_exp,name):
    cur_exp["ftEPR.PlsSPELLISTSlct"].value = name
    
    if cur_exp["ftEPR.PlsSPELLISTSlct"].value != name:
        print("WARNING: Phase cycling did not change")
        return 0
    else:
        return 1
def get_PulseSpel_experiment(cur_exp):
    return cur_exp["ftEPR.PlsSPELEXPSlct"].value
    
def set_PulseSpel_experiment(cur_exp,name):
    cur_exp["ftEPR.PlsSPELEXPSlct"].value = name
    
    if cur_exp["ftEPR.PlsSPELEXPSlct"].value != name:
        print("WARNING: Pulse Spel Experiment did not change")
        return 0
    else:
        return 1

def get_Acquistion_mode(cur_exp):
    return cur_exp["ftEPR.FTAcqModeSlct"].value    
    
def set_Acquistion_mode(cur_exp, mode:int):
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
              
def compile_PulseSpel_prg():
    Xepr.XeprCmds.aqPgShowPrg()
    Xepr.XeprCmds.aqPgCompile()
    time.sleep(0.5)

def compile_PulseSpel_def():
    Xepr.XeprCmds.aqPgShowDef()
    Xepr.XeprCmds.aqPgCompile()
    time.sleep(0.5)


## Section on phase control

class phase():
    """
    A class for the control of phases.
    Options: "cwBridge.SignalPhase", "ftBridge.BrXPhase", "ftBridge.BrYPhase", "ftBridge.BrMinXPhase", "ftBridge.BrMinYPhase"
    """
    def __init__(name:str):
        self.name = name
        self.max = hidden[self.name].aqGetParMaxValue
        self.min = hidden[self.name].aqGetParMinValue
        self.course_step = hidden[self.name].aqGetParCoarseSteps
        self.fine_step = hidden[self.name].aqGetParFineSteps
    
    def get_value(self):
        return hidden[self.name].value
    
    def set_value(self,new_value):
        hidden[self.name.value] = new_value
        return hidden[self.name].value

class attenuator():
    """
    A class for the control of both stepped and variable attenuators
    """
               
               