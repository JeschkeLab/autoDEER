## This is a script for function that run a resonator profile experiment and analyse the data. 

import time
import numpy as np
from autoDeer.hardware import xepr_api
import scipy.fft as fft

def run_nutation(api:xepr_api,pulse_lengths,delays,steps,avgs):
    
    # exp_file = '/home/xuser/Desktop/huka/autoDeer/autoDeer/PulseSpel/autoDEER_4p.exp'
    exp_file = MODULE_DIR + "/PulseSpel/res_pro.exp"
    # def_file = '/home/xuser/Desktop/huka/autoDeer/autoDeer/PulseSpel/autoDEER_4p.def'
    def_file = MODULE_DIR+"/PulseSpel/res_pro.def"

    # 
    api.set_ReplaceMode(False) #Turn replace mode off
    api.set_set_PhaseCycle(True)
    
    # Set the number of points per trace
    d3 = delays[1] - 180
    num_points =  np.floor((delays[1]+delays[2]-d3-200)/steps[0])
    change_DEER_length(exp_file,num_points)
    

    api.set_PulseSpel_exp_filepath(exp_file)
    api.set_PulseSpel_def_filepath(def_file)
    api.compile_PulseSpel_prg()
    api.compile_PulseSpel_def()         
    
    # Set Pulse Lengths
    api.set_PulseSpel_var("p0",pulse_lengths[1])
    api.set_PulseSpel_var("p1",pulse_lengths[0])
    api.set_PulseSpel_var("p2",pulse_lengths[0])
    
    # Set delays
    api.set_PulseSpel_var("d0",delays[0])
    api.set_PulseSpel_var("d1",delays[1])
    api.set_PulseSpel_var("d2",delays[2])
    api.set_PulseSpel_var("d3",d3)
    # Set Steps
    api.set_PulseSpel_var("d30",steps[0])
    api.set_PulseSpel_var("d31",steps[1])
    api.set_PulseSpel_var("d29",steps[2])
    # Set counters
    api.set_PulseSpel_var("h",avgs[0])
    api.set_PulseSpel_var("n",avgs[1])
    api.set_PulseSpel_var("m",avgs[2])
    
    api.set_PulseSpel_experiment("DEER")
    api.set_PulseSpel_phase_cycling("DEER run")
    api.set_Acquistion_mode(1)
    
    # Compile Defs and Program
    api.compile_PulseSpel_prg()
    api.compile_PulseSpel_def()  

    # Run Experiment
    api.run_exp()
    time.sleep(1)
    return 1


def get_nutations(api:xepr_api,channel,nu,field,step,nx:int=128):

    min_freq = nu[0]
    max_freq = nu[1]

    freq_table = np.arange(min_freq,max_freq,step)

    n = len(freq_table)

    if len(field) == 2:
        input_freq = field[0]
        input_field =  field[1]
        start_field = input_field * min_freq / input_freq
    elif len(field) == 1:
        start_field = field
    
    field_table = freq_table * start_field/min_freq
    
    # go to start field /  freq
    api.set_field(field_table[0],hold=True)
    api.set_freq(freq_table[0])

    nut_data = np.zeros((n,nx),dtype=np.complex64)

    for i in range(0,n):
        api.set_field(field_table[i],hold=True)
        api.set_freq(freq_table[i])
        api.run_exp()
        while api.is_exp_running():
            time.sleep(0.5)
        
        dataset = api.acquire_dataset()
        nut_data[i,:] = dataset.data

    return nut_data


def calc_res_prof(nut_data):

    fnut = fft.fft(nut_data)
    fnut = fft.fftshift(fnut)






