from posixpath import expanduser
import numpy as np
import time
import os,sys;
import XeprAPI
import deerlab as dl


# def connect_Xepr():
#     try:
#         Xepr = XeprAPI.Xepr()
#         return Xepr
Xepr = None
cur_exp = None

def set_XEPR_global(Xepr_inst): 
    global Xepr
    Xepr = Xepr_inst

def get_XEPR_global():
    if Xepr != None:
        return Xepr
    else:
        raise RuntimeError("Can't find XEPR instance")


def set_cur_exp_global(cur_exp_inst):
    global cur_exp
    cur_exp = cur_exp_inst

def get_cur_exp_global():
    if cur_exp != None:
        return cur_exp
    else:
        raise RuntimeError("Can't find current experiment")

def is_exp_running(): # Untested
    return Xepr.AQ_EXP_RUNNING()

def acquire_dataset(): # Untested
    """
    This function acquire the dataset
    """
    dataset = cur_exp.XeprDataset()
    size = dataset.size() # This needs to be checked to see what form this precisely comes as
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
    This script detects the end of the scan and acquires the data set. Currently only 1D
    """
    current_scan = cur_exp.getParam("NbScansDone").value
    x_length = int(cur_exp.getParam("XSpecRes").value)
    time_per_point = cur_exp.getParam("ShotRepTime").value *1e-6*cur_exp.getParam("ShotsPLoop").value*2
    trace = np.zeros(x_length, dtype = np.complex64)
    while cur_exp.getParam("NbScansDone").value == current_scan:
        time.sleep(time_per_point)
    time.sleep(time_per_point)
    trace = Xepr.XeprDataset().O
    time_axis = Xepr.XeprDataset().X
    return time_axis,trace
def acquire_scan_at(scan_num:int):
    """
    This script acquires the scan after a spacific number of scans
    """
    x_length = int(cur_exp.getParam("XSpecRes").value)
    time_per_point = cur_exp.getParam("ShotRepTime").value *1e-6*cur_exp.getParam("ShotsPLoop").value*2
    while cur_exp.getParam("NbScansDone").value != scan_num:
        time.sleep(time_per_point*x_length/2)
    return acquire_scan()

  

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
    
def get_PulseSpel_exp_name(cur_exp):
    return os.path.basename(cur_exp["ftEPR.PlsSPELPrgPaF"].value)

def set_PulseSpel_exp_name(cur_exp,fullpath):
    Xepr.XeprCmds.aqPgLoad(fullpath)

def get_PulseSpel_def_name(cur_exp):
    return os.path.basename(cur_exp["ftEPR.PlsSPELGlbPaF"].value) 

def set_PulseSpel_def_name(cur_exp,fullpath):
    Xepr.XeprCmds.aqPgDefLoad(fullpath)

def get_PulseSpel_phase_cycling(cur_exp):
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
              

              