import XeprAPI
import deerlab as dl
import hardware.xepr_api_adv as api
import numpy as np
import os,sys;


def change_DEER_length(path,new_length:int):
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
        
def deerlab_next_scan():
    t,Vexp = api.acquire_scan()
    t = t/1000
    Vexp = dl.correctphase(Vexp)
    t = dl.correctzerotime(Vexp,t)
    r = dl.time2dist(t)          # distance axis, nm
    fit = dl.fitmodel(Vexp,t,r,'P',dl.bg_hom3d,dl.ex_4pdeer,verbose=False)
    sigma = dl.noiselevel(Vexp,'der')
    #print(dl.time2dist(t))
    return [fit, sigma]  

def run_4pDeer(api,pulse_lengths,delays,steps,avgs):
    
    exp_file = '/home/xuser/Desktop/huka/autoDeer/autoDeer/PulseSpel/autoDEER_4p.exp'
    def_file = '/home/xuser/Desktop/huka/autoDeer/autoDeer/PulseSpel/autoDEER_4p.def'
    # 
    api.set_ReplaceMode(False) #Turn replace mode off
    
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

def main_run(api,ps_length:int,delays,filename:str,path:str): # This follows the same structure as the main run in Paramter Optimization
    file = save_file()
    file.create_file(path + filename + ".h5")
    
    run_4pDeer(api,ps_length,delays,steps,avgs)
    
    while api.is_exp_running() == True:
        time.sleep(1)
        
    DEER1_dataset = api.acquire_dataset()
    
    file.save_experimental_data(cp_dataset,"DEER quick")

    api.xepr_save( path+ 'DEER_quick' + filename)
    
    