# This is set of script for running tuning experiments
#
#
# Hugo Karas 2021

import time
import numpy as np
import importlib


MODULE_DIR = importlib.util.find_spec('pydeernet').submodule_search_locations


def get_pulse_trans(api,ps_length,d0):

    # Setting the location of the pulse_spel
    #def_name = '/home/xuser/Desktop/huka/autoDeer/autoDeer/PulseSpel/param_opt.def'
    def_name = MODULE_DIR + '/PulseSpel/param_opt.def'
    #exp_name = '/home/xuser/Desktop/huka/autoDeer/autoDeer/PulseSpel/param_opt.exp'
    exp_name = MODULE_DIR + '/PulseSpel/param_opt.exp'
    
    api.set_ReplaceMode(False) #Turn replace mode off

    # Check that what pulse spel scripts are loaded and compile
    if api.get_PulseSpel_def_filename() != def_name:
        api.set_PulseSpel_def_filepath(def_name)
        time.sleep(0.5)
    if api.get_PulseSpel_exp_filename() != exp_name:
        api.set_PulseSpel_exp_filepath(exp_name)
        time.sleep(0.5)

    api.compile_PulseSpel_prg()
    api.compile_PulseSpel_def() 

    # Set pulse lengths
    api.set_PulseSpel_var("p0",ps_length[0])
    api.set_PulseSpel_var("p1",ps_length[1])

    # Selecting the experiment
    api.set_PulseSpel_experiment("Phase Control")
    api.set_PulseSpel_phase_cycling("None")

    api.set_Acquistion_mode(1) # Run from Pulse Spel

    # Run Experiment
    api.run_exp()
    time.sleep(1)

    return 1

def calc_phase(api):
    d0 = 600
    ps_length = [16,32]

    get_pulse_trans(api,ps_length,d0)

    t,V = api.acquire_dataset()
    time_step = t[1]-t[0]
    V_cor = DC_cor(V)
    real = np.trapz(np.real(V_cor),x=t)
    imag = np.trapz(np.imag(V_cor),x=t)
    phase = real/np.sqrt(real**2 + imag**2)
    return phase

def DC_cor(V):
    num_points = len(V)
    offset =  np.mean(V[int(num_points*0.75):])
    return V - offset
