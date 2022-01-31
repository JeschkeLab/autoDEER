# This is set of script for running tuning experiments
#
#
# Hugo Karas 2021

import time
from turtle import end_fill
from xml.etree.ElementTree import TreeBuilder
import numpy as np
import importlib
from scipy.optimize import minimize_scalar


MODULE_DIR = importlib.util.find_spec('autoDeer').submodule_search_locations


def setup_pulse_trans(api,ps_length:tuple,d0):

    # Setting the location of the pulse_spel
    #def_name = '/home/xuser/Desktop/huka/autoDeer/autoDeer/PulseSpel/param_opt.def'
    def_name = MODULE_DIR[0] + '/PulseSpel/param_opt.def'
    #exp_name = '/home/xuser/Desktop/huka/autoDeer/autoDeer/PulseSpel/param_opt.exp'
    exp_name = MODULE_DIR[0] + '/PulseSpel/param_opt.exp'
    
    api.set_ReplaceMode(True) #Turn replace mode off

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

    return 1

def calc_phase(t,V):
    time_step = t[1]-t[0]
    V_cor = DC_cor(V)
    real = np.trapz(np.real(V_cor),x=t)
    imag = np.trapz(np.imag(V_cor),x=t)
    phase = np.arctan(real/imag)
    # phase = real/np.sqrt(real**2 + imag**2)
    return phase

def DC_cor(V):
    num_points = len(V)
    offset =  np.mean(V[int(num_points*0.75):])
    return V - offset

def tune(api,d0:int,channel:str = 'main',phase_target:str = 'R+'):

    channel_opts = ['main', '+<x>','-<x>','+<y>','-<y>']
    phase_opts = ['R+','R-','I+','I-']
    if not channel in channel_opts:
        raise ValueError(f'Channel must be one of: {channel_opts}')
    
    if not phase_target in phase_opts:
        raise ValueError(f'Phase target must be one of: {phase_opts}')

    setup_pulse_trans(api,[16,32]) #TODO auto-finding of do though the abs pulse

    if channel == 'main':
        lb = 0.0
        ub = 4095.0
        start = 2000.0
        tol = 1
        maxiter = 10
        phase_channel = 'SignalPhase'
    elif channel in ['+<x>','-<x>','+<y>','-<y>']:
        lb = 0.0
        ub = 100.0
        start = 50.0
        tol = 0.1
        maxiter = 10
        if channel == '+<x>':
            phase_channel = 'BrXPhase'
        elif channel == '-<x>':
            phase_channel = 'MinBrXPhase'
        elif channel == '+<y>':
            phase_channel = 'BrYPhase'
        elif channel == '-<y>':
            phase_channel = 'MinBrYPhase'


    if phase_target in ['R+','I+']:
        phase_aim = np.pi / 2
    elif phase_target in ['R-','I-']:
        phase_aim = -np.pi / 2
    
    
    if phase_target in ['I+','I-']:
        imag_target = True
    else:
        imag_target = False
    
    
    def objecive(x,*args):
        '''x is the given phase setting. Args are (phase_target,imag_target)'''
        api.hidden[phase_channel].value = x # Set phase to value
        api.exp_run()
        time.sleep(5)

        t,v = api.acquire_scan()
        v_cor = DC_cor(v)

        if args[1] == True:
            v = -1j * v

        phase = calc_phase(t,v_cor)
        # Calc Phase
        print(f'Phase Setting = {x} - Phase = {phase}')
        return phase - args[0]

    
    output = minimize_scalar(objecive,method='bounded',bounds=[lb,ub],args=(phase_target,imag_target))

    return output.x