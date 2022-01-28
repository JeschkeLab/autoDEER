# This is set of script for running paramter optimization for DEER/PELDOR
# experiments
#
#
# Hugo Karas 2021

# import imp
# from autoDeer.hardware.xepr_api_adv import xepr_api 
import time
import numpy as np
#from scipy.interpolate import RBFInterpolator
from scipy.interpolate import Rbf as RBFInterpolator
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
from autoDeer.TwoD_Experiment import TwoD_Experiment
from autoDeer.File_Saving import save_file
import logging

po_log = logging.getLogger('core.Param_Opt')


def carr_purcell_run(api,ps_length,d0,sweeps=4,steps=100,nuc_mod=[1,1]):
    
    # Setting the location of the pulse_spel
    def_name = '/home/xuser/Desktop/huka/autoDeer/autoDeer/PulseSpel/param_opt.def'
    exp_name = '/home/xuser/Desktop/huka/autoDeer/autoDeer/PulseSpel/param_opt.exp'
    
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
    api.set_PulseSpel_var("p0",ps_length[1])
    api.set_PulseSpel_var("p1",ps_length[0])
    api.set_PulseSpel_var("p2",ps_length[0])

    # Set Pulse Delays
    api.set_PulseSpel_var("d0",d0)
    api.set_PulseSpel_var("d1",200) # Starting pulse delay

    # Set Pulse Steps
    api.set_PulseSpel_var("d12",steps) # Set pulse delay steps
    api.set_PulseSpel_var("d16",nuc_mod[0]) # Nuclear Modulation steps

    # Set Averaging loops
    api.set_PulseSpel_var("h",4) # Shots per points
    api.set_PulseSpel_var("n",sweeps) # Sweeps
    api.set_PulseSpel_var("m",nuc_mod[1]) # Nuclear Modulation sweeps


    # Selecting the experiment
    api.set_PulseSpel_experiment("Carr Purcell exp")
    api.set_PulseSpel_phase_cycling("16_Step")

    api.set_Acquistion_mode(1) # Run from Pulse Spel

    # Run Experiment
    api.run_exp()
    time.sleep(1)

    return 1

def carr_purcell_analysis(dataset):
    """
    This function loads a carr_purcell analysis and finds decay constant
    """
    t = dataset.time
    trace = dataset.data

    trace = np.real(trace)
    trace = trace - np.min(trace)
    trace = trace / np.max(trace)

    rbf = RBFInterpolator(t,trace)
    xi = np.arange(np.min(t),np.max(t),50)
    yi = rbf(xi)

    cp_decay = xi[np.argmin(abs(yi-(1/np.e)))]
    xmax = 2.5 * cp_decay

    # Perform sanity checks on calculated value,
    if xmax < 300:
        print('Carr_purcell_xmax very low: please check!')
        print(f'Carr Purcell Decay estimate is:{xmax}ns')
        print(f'xmax estimate is:{xmax}ns')
        check = None
        while check == None:
            check = input('Type Y to confirm:')
            if check == 'Y' or check == 'y':
                return xmax
            else:
                raise ValueError("xmax is too low")
    elif xmax > 0.8*np.max(t):
        print('Carr_purcell_xmax very high: please check!')
        print(f'xmax estimate is:{xmax}ns')
        check = None
        while check == None:
            check = input('Type Y to confirm:')
            if check == 'Y' or check == 'y':
                return xmax
            else:
                raise ValueError("xmax is too high")
                
    print(f' Carr Purcell Xmax set to {xmax} ns')
        
    return xmax

def carr_purcell_plot(dataset):
    """
    This function plots the carr purcell trace, with 1/e and the calculated max time. 
    """
    t = dataset.time
    trace = dataset.data
    
    fig = plt.figure(figsize=(6,6),dpi=150)
    axs = fig.subplots(1,1)
    
    trace = np.real(trace)
    trace = trace - np.min(trace)
    trace = trace / np.max(trace)

    axs.plot(t,trace,'r')
    axs.plot([np.min(t),np.max(t)],[1/np.e,1/np.e],'k')
    axs.set_xlabel(r'$\tau_1 = \tau_2 / (us)$')
    axs.set_ylabel('Normalized Amplitude')
    axs.set_title('Carr Purcell Experiment')

    return fig

def tau2_scan_run(api,ps_length,d0,tau1):
    
    # Setting the location of the pulse_spel
    def_name = '/home/xuser/Desktop/huka/autoDeer/autoDeer/PulseSpel/param_opt.def'
    exp_name = '/home/xuser/Desktop/huka/autoDeer/autoDeer/PulseSpel/param_opt.exp'
    
    api.set_ReplaceMode(False) #Turn replace mode off
    
    # Check that what pulse spel scripts are loaded and compile
    if api.get_PulseSpel_def_filename() != def_name:
        api.set_PulseSpel_def_filepath(def_name)
    if api.get_PulseSpel_exp_filename() != exp_name:
        api.set_PulseSpel_exp_filepath(exp_name)

    api.compile_PulseSpel_prg()
    api.compile_PulseSpel_def() 
    
    # Set pulse lengths
    api.set_PulseSpel_var("p0",ps_length[1])
    api.set_PulseSpel_var("p1",ps_length[0])
    api.set_PulseSpel_var("p2",ps_length[0])

    # Set Pulse Delays
    api.set_PulseSpel_var("d0",d0)
    api.set_PulseSpel_var("d1",tau1) # Starting pulse delay
    api.set_PulseSpel_var("d2",200) # Starting pulse delay

    # Set Pulse Steps
    api.set_PulseSpel_var("d14",100) # Set pulse delay steps

    # Set Averaging loops
    api.set_PulseSpel_var("h",4) # Shots per points
    api.set_PulseSpel_var("n",4)

    # Selecting the experiment
    api.set_PulseSpel_experiment("tau 2 scan")
    api.set_PulseSpel_phase_cycling("16_Step")

    api.set_Acquistion_mode(1) # Run from Pulse Spel
    
    # Run Experiment
    api.run_exp()
    time.sleep(1)

    return 1

def twoD_scan(api,ps_length,delays,steps,loops):
    # Setting the location of the pulse_spel
    def_name = '/home/xuser/Desktop/huka/autoDeer/autoDeer/PulseSpel/param_opt.def'
    exp_name = '/home/xuser/Desktop/huka/autoDeer/autoDeer/PulseSpel/param_opt.exp'
    
    api.set_ReplaceMode(False) #Turn replace mode off
    
    # Check that what pulse spel scripts are loaded and compile
    if api.get_PulseSpel_def_filename() != def_name:
        api.set_PulseSpel_def_filepath(def_name)
    if api.get_PulseSpel_exp_filename() != exp_name:
        api.set_PulseSpel_exp_filepath(exp_name)

    api.compile_PulseSpel_prg()
    api.compile_PulseSpel_def() 

    # Set pulse lengths
    api.set_PulseSpel_var("p0",ps_length[1])
    api.set_PulseSpel_var("p1",ps_length[0])
    api.set_PulseSpel_var("p2",ps_length[0])

    # Set Pulse Delays
    api.set_PulseSpel_var("d0",delays[0])
    api.set_PulseSpel_var("d1",delays[1]) # Starting pulse delay
    api.set_PulseSpel_var("d2",delays[2]) # Starting pulse delay

    # Set Pulse Steps
    api.set_PulseSpel_var("d12",steps[0]) # tau 1 step
    api.set_PulseSpel_var("d14",steps[1]) # tau 2 step


    # Set Averaging loops
    api.set_PulseSpel_var("h",loops[0]) # Shots per points
    api.set_PulseSpel_var("n",loops[1]) # Sweeps

    # Selecting the experiment
    api.set_PulseSpel_experiment("2D Dec. 64")
    api.set_PulseSpel_phase_cycling("16_Step")

    api.set_Acquistion_mode(1) # Run from Pulse Spel
    
    # Run Experiment
    api.run_exp()
    time.sleep(1)

    
    return 1

def main_run(api,ps_length:int,d0:int,filename:str,path:str):

    file = save_file()
    file.create_file(path + filename + ".h5")
    # Start the carr_purcell_run
    carr_purcell_run(api,ps_length,d0)
    cp_meta = {'Pulse Lengths':ps_length,'d0':d0,'start time':time.strftime("%Y/%m/%d - %H:%M")}
    # Detect when experiments is finished and save data
    while api.is_exp_running() == True:
        time.sleep(1)
    cp_meta.update({'end time':time.strftime("%Y/%m/%d - %H:%M")})
    # Acquire complete data set
        ##cp_t,cp_data = api.acquire_dataset()
    cp_dataset = api.acquire_dataset()
    file.save_experimental_data(cp_dataset,"carr_purcell",meta=cp_meta)

    api.xepr_save( path+ 'cp_q_' + filename)
    # Save complete data set using bruker formats

    # Identify the max time
    cp_max = carr_purcell_analysis(cp_dataset)
    cp_fig = carr_purcell_plot(cp_dataset)
    cp_fig.show()
    po_log.info(f'Carr-Purcell decay maximum of {cp_max}ns')

    # Now that the cp max time has been caclulcated we need to repeat this along
    # the vertical axis. This might be unesessary, but they will likely be differnt.

    # Run the tau2_scan
    tau1 = 400 #ns
    tau2_scan_run(api,ps_length,d0,tau1)
    tau2_meta = {'Pulse Lengths':ps_length,'d0':d0,'tau1':tau1,'start time':time.strftime("%Y/%m/%d - %H:%M")}

    # Wait until tau2 scan finishes
    while api.is_exp_running() == True:
        time.sleep(1)
    tau2_meta.update({'end time':time.strftime("%Y/%m/%d - %H:%M")})

    # Acquire complete data set
    # tau2_t,tau2_data = api.acquire_dataset()
    tau2_dataset = api.acquire_dataset()

    file.save_experimental_data(tau2_dataset,"tau2_scan",meta=tau2_meta)
    api.xepr_save( path+ 'tau2_' + filename)

    
    # Identify the max time
    tau2_max = carr_purcell_analysis(tau2_dataset)
    tau2_fig = carr_purcell_plot(tau2_dataset)
    tau2_fig.show()
    po_log.info(f'Tau2 decay maximum of {tau2_max}ns')

    # Now that the two maxes have been found take the largest one of these to be the max trace.
    max_tau = max([cp_max,tau2_max])
    print(f" Maximum tau2 set to {int(max_tau)} ns")
    po_log.info(f'Overall decay maximum of {max_tau}ns')
    
    tau_step = np.floor((max_tau - 200)/64)
    tau_step = (tau_step//2)*2 # make sure this is a multiple of 2 for c-floor
    delays = [d0,200,200]
    steps = [tau_step,tau_step]
    loops = [4,4]
    twoD_scan(api,ps_length,delays,steps,loops)
    twoD_meta = {'Pulse Lengths':ps_length,'delays':delays,'steps':steps,'loops':loops,'start time':time.strftime("%Y/%m/%d - %H:%M")}

    
    while api.is_exp_running() == True:
        time.sleep(1)
    twoD_meta.update({'end time':time.strftime("%Y/%m/%d - %H:%M")})



    twoD_dataset = api.acquire_scan_2d()
    file.save_experimental_data(twoD_dataset,"2D_exp",meta=twoD_meta)
    api.xepr_save( path+ '2D_dec_' + filename)

    last_scan = TwoD_Experiment()
    last_scan.import_dataset(twoD_dataset)
    last_scan.snr_normalize(loops[0]*loops[1]*16)
    last_scan.calculate_optimal()
    print(f'The optimal pulse delays for 4p DEER are: {last_scan.time_4p}')
    po_log.info(f'Optimal pulse delays for 4p DEER are: {last_scan.time_4p} ns')

    # TODO: ADD uncertainty estimation to 2D plot
    # Once more than one scan has been collected look into how the uncertainty can be estimated
    # This is significantly harder 
    
    
    return 1
