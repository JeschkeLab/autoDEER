# This is set of script for running paramter optimization for DEER/PELDOR
# experiments
#
#
# Hugo Karas 2021

import xepr_api_adv as api
import time
import numpy as np
from scipy.interpolate import RBFInterpolator
def carr_purcell_run(cur_exp,ps_length,d0):
    
    # Setting the location of the pulse_spel
    def_name = 'param_opt.def'
    exp_name = 'param_opt.exp'
    
    api.set_ReplaceMode(cur_exp,False) #Turn replace mode off
    
    # Check that what pulse spel scripts are loaded and compile
    if api.get_PulseSpel_def_name() != def_name:
        api.set_PulseSpel_def_name(cur_exp,def_name)
    if api.get_PulseSpel_exp_name() != exp_name:
        api.set_PulseSpel_def_name(cur_exp,exp_name)

    api.compile_PulseSpel_prg()
    api.compile_PulseSpel_def() 

    # Set pulse lengths
    api.set_PulseSpel_var(cur_exp,"p0",ps_length)
    api.set_PulseSpel_var(cur_exp,"p1",ps_length)
    api.set_PulseSpel_var(cur_exp,"p2",ps_length)

    # Set Pulse Delays
    api.set_PulseSpel_var(cur_exp,"d0",d0)
    api.set_PulseSpel_var(cur_exp,"d1",200) # Starting pulse delay

    # Set Pulse Steps
    api.set_PulseSpel_var(cur_exp,"d30",100) # Set pulse delay steps

    # Set Averaging loops
    api.set_PulseSpel_var(cur_exp,"h",4) # Shots per points
    api.set_PulseSpel_var(cur_exp,"n",2) # Sweeps

    # Selecting the experiment
    api.set_PulseSpel_experiment(cur_exp,"carr_purcell")
    api.set_PulseSpel_phase_cycling(cur_exp,"16_Step")

    api.set_Acquistion_mode(cur_exp,1) # Run from Pulse Spel

    # Run Experiment
    cur_exp.aqExpRun()

    return 1

def carr_purcell_analysis(time,trace):
    """
    This function loads a carr_purcell analysis and finds decay constant
    """
    trace = np.real(trace)
    trace = trace - np.min(trace)
    trace = trace / np.max(trace)

    rbf = RBFInterpolator(time,trace)
    xi = np.arange(np.min(time),np.max(time),50)
    yi = rbf(xi)

    cp_decay = xi[np.argmin(abs(yi-(1/np.e)))]
    xmax = 6 * cp_decay

    # Perform sanity checks on calculated value,
    if xmax < 300:
        print('Carr_purcell_xmax very low: please check!')
        print(f'Carr Purcell Decay estimate is:{xmax}ns')
        print(f'xmax estimate is:{xmax}ns')
        check = input('Type Y to confirm:')
        if check == 'Y' or check == 'y':
            return xmax
        else:
            raise ValueError("xmax is too low")
    elif xmax > 0.8*np.max(time):
        print('Carr_purcell_xmax very high: please check!')
        print(f'xmax estimate is:{xmax}ns')
        check = input('Type Y to confirm:')
        if check == 'Y' or check == 'y':
            return xmax
        else:
            raise ValueError("xmax is too high")
        
    return xmax

def carr_purcell_plot(time,trace):
    """
    This function plots the carr purcell trace, with 1/e and the calculated max time. 
    """
    return 0 

def tau2_scan(cur_exp,tau1):
    api.set_ReplaceMode(cur_exp,False) #Turn replace mode off
    return 0

def twoD_scan(cur_exp):
    api.set_ReplaceMode(cur_exp,False) #Turn replace mode off
    return 0

def main_run(ps_length,d0):

    cur_exp = api.get_cur_exp_global()

    # Start the carr_purcell_run
    carr_purcell_run(cur_exp,ps_length,d0)

    # Detect when experiments is finished and save data
    while api.is_exp_running() == True:
        time.sleep(1)
    
    # Acquire complete data set
    cp_t,cp_data = api.acquire_dataset()
    # Save complete data set using bruker formats

    # Identify the max time


    # Run the tau2_scan
    while api.is_exp_running() == True:
        time.sleep(1)
    return 0

