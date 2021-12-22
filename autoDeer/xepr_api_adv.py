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

def deerlab_next_scan():
    t,Vexp = acquire_scan()
    t = t/1000
    Vexp = dl.correctphase(Vexp)
    t = dl.correctzerotime(Vexp,t)
    r = dl.time2dist(t)          # distance axis, nm
    fit = dl.fitmodel(Vexp,t,r,'P',dl.bg_hom3d,dl.ex_4pdeer,verbose=False)
    sigma = dl.noiselevel(Vexp,'der')
    #print(dl.time2dist(t))
    return [fit, sigma]    

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
              
# These are ETH DEER Specific function. Will be moved to a seperate file at a later data
def change_DEER_length(cur_exp,path,new_length:int):
    """This is a HIGHLY specific function to change  a line in a specific pulse spel file. This will break your pulse spel script if applied to any other file."""
    with open(path, 'r') as file:
        data = file.readlines()
    if "dim4" in data[16]:
        data[16] = f' dim4 s[{int(new_length)}]              ;     dimension [sx] for DEER\n'
        print(f'DEER Dimension changed to {int(new_length)}')
        with open(path, 'w') as file:
            file.writelines( data )
    else:
        print("ERROR: Can't update the DEER Length. Check pulse Spel Experiment file")
        

def run_4pDeer(cur_exp,pulse_lengths,delays,steps,avgs):
    
    exp_file = '/home/xuser/Desktop/huka/autoDeer/autoDeer/PulseSpel/autoDEER_4p.exp'
    def_file = '/home/xuser/Desktop/huka/autoDeer/autoDeer/PulseSpel/autoDEER_4p.def'
    # 
    set_ReplaceMode(cur_exp,False) #Turn replace mode off
    
    # Set the number of points per trace
    d3 = delays[1] - 180
    num_points =  np.floor((delays[1]+delays[2]-d3-200)/steps[0])
    change_DEER_length(cur_exp,exp_file,num_points)
    

    set_PulseSpel_exp_name(cur_exp,exp_file)
    set_PulseSpel_def_name(cur_exp,def_file)
    compile_PulseSpel_prg()
    compile_PulseSpel_def()         
    
    # Set Pulse Lengths
    set_PulseSpel_var(cur_exp,"p0",pulse_lengths)
    set_PulseSpel_var(cur_exp,"p1",pulse_lengths)
    set_PulseSpel_var(cur_exp,"p2",pulse_lengths)
    
    # Set delays
    set_PulseSpel_var(cur_exp,"d0",delays[0])
    set_PulseSpel_var(cur_exp,"d1",delays[1])
    set_PulseSpel_var(cur_exp,"d2",delays[2])
    set_PulseSpel_var(cur_exp,"d3",d3)
    # Set Steps
    set_PulseSpel_var(cur_exp,"d30",steps[0])
    set_PulseSpel_var(cur_exp,"d31",steps[1])
    set_PulseSpel_var(cur_exp,"d29",steps[2])
    # Set counters
    set_PulseSpel_var(cur_exp,"h",avgs[0])
    set_PulseSpel_var(cur_exp,"n",avgs[1])
    set_PulseSpel_var(cur_exp,"m",avgs[2])
    
    set_PulseSpel_experiment(cur_exp,"DEER")
    set_PulseSpel_phase_cycling(cur_exp,"DEER run")
    set_Acquistion_mode(cur_exp,1)
    
    # Compile Defs and Program
    compile_PulseSpel_prg()
    compile_PulseSpel_def()  

    # Run Experiment
    cur_exp.aqExpRun()
    return 1
              