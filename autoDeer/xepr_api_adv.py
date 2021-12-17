import numpy as np
import time
import os,sys;
import XeprAPI


def acquire_scan():
    """
    This script detects the end of the scan and acquires the data set. Currently only 1D
    """
    current_scan = cur_exp.getParam("NbScansDone").value
    time_per_point = cur_exp.getParam("ShotRepTime").value *1e-6*cur_exp.getParam("ShotsPLoop").value*2
    trace = np.zeros(x_length, dtype = np.complex64)
    while cur_exp.getParam("NbScansDone").value == current_scan:
        time.sleep(time_per_point)
    time.sleep(time_per_point)
    trace = Xepr.XeprDataset().O
    trace = Xepr.XeprDataset().X
    return time,trace
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
    Vexp,t = acquire_scan()
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
    cur_exp["ftEXP.PlsSPELSetVar"].value = str(variable) + " = " str(int(value))
    
def set_ReplaceMode(cur_exp,state=False)
    if state:
        value = 'On'
        print('DANGER: Replace Mode turned ON!')
    else:
        value = 'Off'
    cur_exp["ftEXP.ReplaceMode"].value = value
    
def get_PulseSpel_exp_name(cur_exp):
    return os.path.basename(cur_exp["ftEPR.PlsSPELPrgPaF"].value)

def get_PulseSpel_def_name(cur_exp):
    return os.path.basename(cur_exp["ftEPR.PlsSPELGlb"].value) 