# This is set of script for running tuning experiments
#
#
# Hugo Karas 2021

import time
import numpy as np
import importlib
from scipy.optimize import minimize_scalar


MODULE_DIR = importlib.util.find_spec('autoDeer').submodule_search_locations


def setup_pulse_trans(api,ps_length:tuple,d0,channel:str = "MainPhase"):

    # Setting the location of the pulse_spel
    #def_name = '/home/xuser/Desktop/huka/autoDeer/autoDeer/PulseSpel/param_opt.def'
    def_name = MODULE_DIR[0] + '/PulseSpel/phase_set.def'
    #exp_name = '/home/xuser/Desktop/huka/autoDeer/autoDeer/PulseSpel/param_opt.exp'
    exp_name = MODULE_DIR[0] + '/PulseSpel/phase_set.exp'
    
    api.set_ReplaceMode(False) #Turn replace mode off
    api.set_PhaseCycle(False)
    
    # Set Specjet No. of Points
    
    api.hidden['NoOfPoints'].value = 512
    api.hidden['NoOfAverages'].value = 20
    api.hidden['NOnBoardAvgs'].value = 60

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
    
    api.set_PulseSpel_var("d0",d0)
    api.set_PulseSpel_var("d1",400)
    api.set_PulseSpel_var("d2",400) # Not required but done anyway

    # Selecting the experiment
    api.set_PulseSpel_experiment("Hahn Echo")
    api.set_PulseSpel_phase_cycling(channel)
    
    # Load setting to tables
    api.set_Acquistion_mode(1)
    api.run_exp()
    time.sleep(5)
    api.stop_exp()
    
    # Start Transient Mode
    api.set_Acquistion_mode(3)
    api.run_exp()
    
    time.sleep(5)
    
    # Read Transient Monde
    api.set_Acquistion_mode(2) # Run from om Pulse Spel usin
    api.run_exp()
    return 1

def calc_phase(t,V):
    time_step = t[1]-t[0]
    V_cor = DC_cor(V)
    real = np.trapz(np.real(V_cor),x=t)
    imag = np.trapz(np.imag(V_cor),x=t)
    phase = np.arctan2(imag,real)
    # phase = real/np.sqrt(real**2 + imag**2)
    return phase

def DC_cor(V):
    num_points = len(V)
    offset =  np.mean(V[int(num_points*0.75):])
    return V - offset

def tune(api,d0:int = 600,channel:str = 'main',phase_target:str = 'R+',ps_length:int = 32):

    if ps_length%4 != 0:
        raise ValueError("ps_length must be a multiple of 4")
    elif ps_length < 4:
        raise ValueError("ps_length must be a greater than or equal to 4")
    else:
        ps_lengths = (int(ps_length/2),int(ps_length))
    
    channel_opts = ['main', '+<x>','-<x>','+<y>','-<y>']
    phase_opts = ['R+','R-','I+','I-']
    if not channel in channel_opts:
        raise ValueError(f'Channel must be one of: {channel_opts}')
    
    if not phase_target in phase_opts:
        raise ValueError(f'Phase target must be one of: {phase_opts}')


    if channel == 'main':
        lb = 0.0
        ub = 4095.0
        start = 2000.0
        tol = 5
        maxiter = 10
        phase_channel = 'SignalPhase'
    elif channel in ['+<x>','-<x>','+<y>','-<y>']:
        lb = 0.0
        ub = 100.0
        start = 50.0
        tol = 1
        maxiter = 10
        if channel == '+<x>':
            phase_channel = 'BrXPhase'
        elif channel == '-<x>':
            phase_channel = 'BrMinXPhase'
        elif channel == '+<y>':
            phase_channel = 'BrYPhase'
        elif channel == '-<y>':
            phase_channel = 'BrMinYPhase'

    setup_pulse_trans(api,(16,32),d0,channel=phase_channel) #TODO auto-finding of do though the abs pulse

    # if phase_target in ['R+','I+']:
    #     phase_aim = np.pi / 2
    # elif phase_target in ['R-','I-']:
    #     phase_aim = -np.pi / 2
    
    if phase_target == 'R+':
        phase_aim = 0
        neg_aim = False
    elif phase_target == 'I+':
        phase_aim = np.pi/2
        neg_aim = False
    elif phase_target == 'I-':
        phase_aim = -np.pi/2
        neg_aim = False

    if phase_target == 'R-':
        phase_aim = 0
        neg_aim = True
    
    
    # if phase_target in ['I+','I-']:
    #     imag_aim = True
    # else:
    #     imag_aim = False
    
    
    def objecive(x,*args):
        '''x is the given phase setting. Args are (phase_target,imag_target)'''
        api.hidden[phase_channel].value = x # Set phase to value
        time.sleep(2)
        api.run_exp()
        while api.is_exp_running():
            time.sleep(1)

        data = api.acquire_scan()
        t = data.time
        v = data.data
        v_cor = DC_cor(v)

        if args[1] == True:
            v_cor = -1 * v_cor

        phase = calc_phase(t,v_cor)
        # Calc Phase
        print(f'Phase Setting = {x:.1f} \t Phase = {phase:.2f} \t Phase Dif = {np.abs(phase - args[0]):.2f}')
        return np.abs(phase - args[0])

    print(f'Phase Aim = {phase_aim:.3f}')
    api.hidden['AverageStart'].value = True

    output = minimize_scalar(objecive,method='bounded',bounds=[lb,ub],args=(phase_aim,neg_aim),options={'xatol':tol,'maxiter':30})
    
    print(f"Optimal Phase Setting for {phase_channel} is: {output.x:.1f}")
    api.hidden[phase_channel].value = output.x
    
    return output.x



def tune_power(api,d0:int = 600,channel:str = 'main',ps_target:int = 32) -> float:

    
    if ps_target%4 != 0:
        raise ValueError("ps_target must be a multiple of 4")
    elif ps_target < 4:
        raise ValueError("ps_target must be a greater than or equal to 4")
    else:
        ps_lengths = (int(ps_target/2),int(ps_target))
    
    
    channel_opts = ['main', '+<x>','-<x>','+<y>','-<y>']
    if not channel in channel_opts:
        raise ValueError(f'Channel must be one of: {channel_opts}')



    if channel == 'main':
        lb = 0.0
        ub = 60.0
        tol = 1
        atten_channel = 'PowerAtten'
    elif channel in ['+<x>','-<x>','+<y>','-<y>']:
        lb = 0.0
        ub = 100.0
        tol  = 2
        if channel == '+<x>':
            atten_channel = 'BrXAmp'
            phase_channel = 'BrXPhase'
        elif channel == '-<x>':
            atten_channel = 'BrMinXAmp'
            phase_channel = 'BrMinXPhase'            
        elif channel == '+<y>':
            atten_channel = 'BrYAmp'
            phase_channel = 'BrYPhase'
        elif channel == '-<y>':
            atten_channel = 'BrMinYAmp'
            phase_channel = 'BrMinYPhase'

            
    setup_pulse_trans(api,ps_lengths,d0,channel = phase_channel) #TODO auto-finding of do though the abs pulse

    def objecive(x,*args):
        '''x is the given attenuator setting'''
        api.hidden[atten_channel].value = x # Set phase to value
        time.sleep(2)
        api.run_exp()
        while api.is_exp_running():
            time.sleep(1)

        data = api.acquire_scan()
        t = data.time
        v = data.data
        v_cor = DC_cor(v)

        inte = -1 *  np.trapz(np.abs(np.real(v_cor)) + np.abs(np.imag(v_cor))) 

        # Calc Phase
        print(f'Attenuator Setting = {x:.1f} \t inte = {inte:.2f} ')
        return inte
    
    api.hidden['AverageStart'].value = True
    output = minimize_scalar(objecive,method='bounded',bounds=[lb,ub],options={'xatol':tol,'maxiter':30})

    print(f"Optimal Attenuator Setting for {atten_channel} is: {output.x:.1f}")
    api.hidden[atten_channel].value = output.x
    return output.x
    


def calc_d0(api,acq_length):
    guess = 500
    
    
    data = xepr.acquire_dataset()
    t = data.time
    d = data.data
    d0 = guess + round(t[np.argmax(abs(d))]/2)*2 - acq_length/2
    return d0