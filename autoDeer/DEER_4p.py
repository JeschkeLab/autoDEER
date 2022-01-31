import importlib
import deerlab as dl
import numpy as np
import os,sys;

from autoDeer.File_Saving import save_file
import time
import importlib

MODULE_DIR = importlib.util.find_spec('autoDeer').submodule_search_locations

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
        
def deerlab_next_scan(api):
    t,Vexp = api.acquire_scan()
    t = t/1000
    Vexp = dl.correctphase(Vexp)
    t = dl.correctzerotime(Vexp,t)
    r = dl.time2dist(t)          # distance axis, nm
    fit = dl.fitmodel(Vexp,t,r,'P',dl.bg_hom3d,dl.ex_4pdeer,verbose=False)
    sigma = dl.noiselevel(Vexp,'der')
    #print(dl.time2dist(t))
    return [fit, sigma]  

def std_deerlab(t,Vexp):
    t = t/1000
    Vexp = dl.correctphase(Vexp)
    t = dl.correctzerotime(Vexp,t)
    r = dl.time2dist(t)          # distance axis, nm
    fit = dl.fitmodel(Vexp,t,r,'P',dl.bg_hom3d,dl.ex_4pdeer,verbose=False)
    sigma = dl.noiselevel(Vexp,'der')
    #print(dl.time2dist(t))
    return [fit, sigma]  

def run_4pDeer(api,pulse_lengths,delays,steps,avgs):
    
    # exp_file = '/home/xuser/Desktop/huka/autoDeer/autoDeer/PulseSpel/autoDEER_4p.exp'
    exp_file = MODULE_DIR + "/PulseSpel/autoDEER_4p.exp"
    # def_file = '/home/xuser/Desktop/huka/autoDeer/autoDeer/PulseSpel/autoDEER_4p.def'
    def_file = MODULE_DIR+"/PulseSpel/autoDEER_4p.def"

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

def main_run(api,ps_length:int,delays,filename:str,path:str,steps = [12,2,2], exp_time=2): # This follows the same structure as the main run in Paramter Optimization
    file = save_file()
    file.open_file(path + filename + ".h5")

    meta_4pDeer = {'Pulse Lengths':ps_length,'d0':delays,'steps':steps,'start time':time.strftime("%Y/%m/%d - %H:%M")}
    run_4pDeer(api,ps_length,delays,steps,[10,2000,1])
    
    time.sleep(exp_time*60*60) # This implements the time limit
    api.stop_exp()

    while api.is_exp_running() == True:
        time.sleep(1)
    meta_4pDeer.update({'end time':time.strftime("%Y/%m/%d - %H:%M")})


    DEER1_dataset = api.acquire_dataset()
    
    dset_raw = file.save_experimental_data(DEER1_dataset,"DEER quick",meta=meta_4pDeer)

    api.xepr_save( path+ 'DEER_quick' + filename)
    
    [fit, sigma] = std_deerlab(DEER1_dataset.time,DEER1_dataset.data)
    dl_meta = {'sigma':sigma,'ddparam':fit.ddparam,'bgparam':fit.bgparam,'exparam':fit.exparam}
    dset_deer = file.save_experimental_data(DEER1_dataset,"DEER quick",dset_name="Distance Distribution",meta=dl_meta)

    return fit, sigma