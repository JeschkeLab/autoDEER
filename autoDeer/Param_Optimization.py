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
import importlib

po_log = logging.getLogger('core.Param_Opt')

MODULE_DIR = importlib.util.find_spec('autoDeer').submodule_search_locations


def carr_purcell_run(api,ps_length,d0,sweeps=4,steps=100,nuc_mod=[1,1]):
    
    # Setting the location of the pulse_spel
    def_name = '/home/xuser/Desktop/huka/autoDeer/autoDeer/PulseSpel/param_opt.def'
    def_name = MODULE_DIR + '/PulseSpel/param_opt.def'
    exp_name = '/home/xuser/Desktop/huka/autoDeer/autoDeer/PulseSpel/param_opt.exp'
    exp_name = MODULE_DIR + '/PulseSpel/param_opt.exp'
    
    api.set_ReplaceMode(False) #Turn replace mode off
    api.set_set_PhaseCycle(True)

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
    def_name = MODULE_DIR + '/PulseSpel/param_opt.def'
    exp_name = '/home/xuser/Desktop/huka/autoDeer/autoDeer/PulseSpel/param_opt.exp'
    exp_name = MODULE_DIR + '/PulseSpel/param_opt.exp'
    
    api.set_ReplaceMode(False) #Turn replace mode off
    api.set_set_PhaseCycle(True)

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
    def_name = MODULE_DIR + '/PulseSpel/param_opt.def'
    exp_name = '/home/xuser/Desktop/huka/autoDeer/autoDeer/PulseSpel/param_opt.exp'
    exp_name = MODULE_DIR + '/PulseSpel/param_opt.exp'
    
    api.set_ReplaceMode(False) #Turn replace mode off
    api.set_set_PhaseCycle(True)
    
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



def run_EDFS(api, d0:int, filename:str, path:str, sweep_size:int = 300, scans:int = 4, shots:int = 10) -> None:
    """run_EDFS Runs an Echo Detected Field Sweep (EDFS)

    Parameters
    ----------
    api : [type]
        The class that talks to the spectrometer this is being applied to.
    d0 : int
        The spectrometer offset delay, d0
    filename : str
        The name of the save file.
    path : str
        The *full* path to the save file
    sweep_size : int, optional
        The number of data points per sweep, by default 300
    scans : int, optional
        The number of scans to be done, by default 4
    shots : int, optional
        The number of shots per point, by default 10
    """
    meta = {'Experiment':'EDFS','D0':d0,'scans':scans,'shots':shots}

    # Open/Create File

    file = save_file()
    file.create_file(path + filename + ".h5")

    # Get Current Frequency
    freq = api.get_freq()

    # Get Current Magnet Value
    field = api.get_field()

    sweep_width = api.set_sweep_width(sweep_size)

    meta.update({'MW Freq':freq,'Field':field,'Sweep Width':sweep_width})
    # Run pulse Spel EDFS

    run_ps_script(api,'/PulseSpel/param_opt',"EDFS","EDFS",[16,32],[d0,400,400],[0,0],[shots,scans])

    while api.is_exp_running() == True:
        time.sleep(1)

    EDFS_data = api.acquire_dataset()
    api.xepr_save(path+ 'EDFS_' + filename)

    # Find the main peak
    peak = EDFS_data.time[np.argmax(EDFS_data.data)]
    
    meta.update({'EDFS_peak':peak})

    file.save_experimental_data(EDFS_data,"EDFS",meta=meta)


    pass

def run_ps_script(api,ps_script:str,exp:str,phase_cycle:str,ps_length:list[int,int],
        delays:list[int,int,int],steps:list[int,int],loops:list[int,int], shrt:int=6000,**kwargs) -> None:
    """run_ps_script Runs a general pulse spell script according to the set paramters. The script must conform to the standard setup.

    Parameters
    ----------
    api : [type]
        The class that talks to the spectrometer this is being applied to.
    ps_script : str
        The *local* path to the pulse spel script
    exp : str
        The name of the experiment to be run.
    phase_cycle : str
        The name of the phase cycle to be run.
    ps_length : list[int,int]
        The length (in ns) of the pulses used [pi/2,pi]
    delays : list[int,int,int]
        The length (in ns) of the delays [d0,d1,d2].
    steps : list[int,int]
        The length (in ns) of the steps sizes [d11,d21]/
    loops : list[int,int]
        The number of iterations in each loop [shots per point,sweeps]
    shrt : int
        The shot repetition time in us, by default 6000us.
    """
    def_name = MODULE_DIR + ps_script + '.def'
    exp_name = MODULE_DIR + ps_script + '.exp' 

    if "ReplaceMode" in kwargs:
        api.set_ReplaceMode(kwargs["ReplaceMode"])
    else:
        api.set_ReplaceMode(False) #Turn replace mode off
    
    if "PhaseCycle" in kwargs:
        api.set_PhaseCycle(kwargs["PhaseCycle"]) 
    else:
        api.set_PhaseCycle(True) #This is the default value

    if api.get_PulseSpel_def_filename() != def_name:
        api.set_PulseSpel_def_filepath(def_name)
    if api.get_PulseSpel_exp_filename() != exp_name:
        api.set_PulseSpel_exp_filepath(exp_name)

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


    api.compile_PulseSpel_prg()
    api.compile_PulseSpel_def() 

    # Selecting experiments
    api.set_PulseSpel_experiment(exp)
    api.set_PulseSpel_phase_cycling(phase_cycle)

    api.set_Acquistion_mode(1) # Run from Pulse Spel

    api.run_exp()
    time.sleep(1)

    pass


# Function to autophase the MPFU channels

class MaxIterError(RuntimeError):
    pass

class BoundError(RuntimeError):
    pass

def secant_method(fun:function,x0:float,x1:float) -> float:
    """ The secant method """
    x2 = (x0*fun(x1) - x1*fun(x0))/(fun(x1)-fun(x0))
    return x2
    
def root_finder(obj:function,x0:float,x1:float,bounds:list[float,float]=[0,1],MaxIter:int=20,tol:float=1,MaxStep:float =None) -> float:
    """root_finder Optimizes the objective function through a modified secant method. 
    This introduces bounds on the method, as well as a maximum number of iterations and a maximum step size. 

    Parameters
    ----------
    obj : function
        The objective function for which the root is to be found. This should only have one input variable, and the root should be located at y=0.
    x0 : float
        The lower initial guess
    x1 : float
        The upper initial guess
    bounds : list[float,float], optional
        The lower and upper bound of the stepper motor. If the optimization leaves these limits it will return a BoundError, by default [0,1]
    MaxIter : int, optional
        The maximum number of steps/iterations per cycle, by default 30
    tol : float, optional
        The condition for when to stop. When the difference between x1 and x0 is below this value it will stop, by default 1
    MaxStep : float, optional
        The maximum step size per step of the secant method, by default None
    Returns
    -------
    float
        The root of this objective function under these conditions

    Raises
    ------
    BoundError
        The secant method moved outside the bounds. Try using a smaller MaxStep size.
    MaxIterError
        The optimization stopped due too many iterations. Try increasing the tolerance or the MaxIter. 
    """
    if x0 == x1:
        x0 = x0-1
        x1 = x1+1
    count = 0
    Found = False
    while (not Found) and (count < MaxIter):
        print(x0,x1)
        x2 = obj(x0,x1)
        
        if (MaxStep != None):
            if (np.abs(x2-x1) > MaxStep):
                sign = np.sign(x2-x1)
                x2 = x1 + sign * MaxStep
        
        if x2 < bounds[0] or x2 > bounds[1]:
            raise BoundError('Outside limits. Try reducing the Maximum Step Size')
        
        x0 = x1
        x1 = x2
        count = count + 1
        diff = np.abs(x1-x0)
        if diff < tol:
            Found = True

    if count == MaxIter:
        print('Warning: Max Iteration reached. Increase tolerence or averaging')
        raise MaxIterError('MaxIter Reached')
    else:
        return (x0+x1)/2

def trans_angle(api,x,channel:str)->float:
    """trans_angle Sets the phase of the appropraite channel, records a dataset, takes the mean and then returns the complex phase.

    Parameters
    ----------
    api : [type]
        The class that talks to the spectrometer this is being applied to.
    x : int or float
        The setting for the phase stepper motor
    channel : str
        Name of channel to be tunned, using XeprAPI notations. Options:['SignalPhase','BrXPhase',MinBrXPhase','BrYPhase','MinBrYPhase']

    Returns
    -------
    float
        The current phase of this chosen channel.
    """
    api.hidden[channel].value = x
    time.sleep(10)
    dataset = api.acquire_dataset()
    val = np.mean(dataset.data)
    angle = np.arctan2(np.imag(val),np.real(val))
    return angle

def tune(api,target:str,channel:str,lim:list,tol:float,MaxStep:float,MaxIter:int=30) -> float:
    """tune Optimizes the phase of a given channel to a specific target position. 
    In the case of the optimization failing due to too many iterations the tolerance is automatically doubled. This is done a maximum 
    of two times.

    Parameters
    ----------
    api : [type]
        The class that talks to the spectrometer this is being applied to.
    target : str
        How to tune the channel. Options:['R+','R-','I+','I-']
    channel : str
        Name of channel to be tunned, using XeprAPI notations. Options:['SignalPhase','BrXPhase',MinBrXPhase','BrYPhase','MinBrYPhase']
    lim : list[float,float]
        The lower and upper bound of the stepper motor. If the optimization leaves these limits it will return a BoundError, by default [0,1]
    tol : float
        The condition for when to stop. When the difference between x1 and x0 is below this value it will stop, by default 1
    MaxStep : float
        The maximum step size per step of the secant method.
    MaxIter : int, optional
        The maximum number of steps/iterations per cycle, by default 30

    Returns
    -------
    float
        The optimal setting for that channel

    Raises
    ------
    ValueError
        If the target isn't one of the options.
    """
    if target == 'I+':
        gaus_int_imag = lambda x: trans_angle(api,x,channel) - np.pi/2
    elif target == 'I-':
        gaus_int_imag = lambda x: trans_angle(api,x,channel) + np.pi/2
    elif target == 'R+':
        gaus_int_imag = lambda x: trans_angle(api,x,channel)
    elif target == 'R-':
        gaus_int_imag = lambda x: trans_angle(api,x,channel) + np.pi
    else:
        raise ValueError("Target must be one of:['R+','R-','I+','I-']")

    if channel == 'SignalPhase':
        obj = lambda x0,x1: round(secant_method(gaus_int_imag,x0,x1))
    else:
        obj = lambda x0,x1: round(secant_method(gaus_int_imag,x0,x1),2)


    try:
        x = root_finder(obj,0.45*lim[1],0.55*lim[1],bounds=lim,tol=tol,MaxIter=MaxIter,MaxStep=MaxStep)
    except MaxIterError:
        tol = tol * 2 
        print(f'Tolerance increased to {tol}')
        try:
            x = root_finder(obj,0.45*lim[1],0.55*lim[1],bounds=lim,tol=tol,MaxIter=MaxIter,MaxStep=MaxStep)
        except MaxIterError:
            tol = tol  * 2
            print(f'Tolerance increased to {tol}') 
            x = root_finder(obj,0.45*lim[1],0.55*lim[1],bounds=lim,tol=tol,MaxIter=MaxIter,MaxStep=MaxStep)
    return x

def MpfuTune(api,target:str,channel:str) -> float:
    """MpfuTune: Tunes the given MPFU channel, such that is closely matches the target condition. This is done through a modified secant method.
    
    
    Parameters
    ----------
    api : [type]
        The class that talks to the spectrometer this is being applied to.
    target : str
         How to tune the channel. Options:['R+','R-','I+','I-']
    channel : str
        Name of channel to be tunned, using XeprAPI notations. Options:['+<X>','-<X>','+<Y>','-<Y>']

    Returns
    -------
    float
        The optimal setting for that channel
    """
    if channel == '+<X>':
        channel = 'BrXPhase'
    elif channel == '-<X>':
        channel == 'MinBrXPhase'
    elif channel == '+<Y>':
        channel == 'BrYPhase'
    elif channel == '-<Y>':
        channel == 'MinBrYPhase'
    
    return tune(api,target,channel,[0,100],0.5,20)

def MainTune(api,target:str) -> int:
    """MainTune  Tunes the global main channel, such that is closely matches the target condition. 
    This is done through a modified secant method, under integer conditions.

    Parameters
    ----------
    api : [type]
        [description]
    target : str
         How to tune the channel. Options:['R+','R-','I+','I-']

    Returns
    -------
    int
        The optimal setting for that channel
    """
    return round(tune(api,target,'SignalPhase',[0,4095],5,100))

def runTune(api,channel:str,d0:int):
    run_ps_script(api,'/PulseSpel/phase_set','Hahn Echo',channel,[16,32],[d0,400,400],[0,0,0],[10,2000,0],ReplaceMode=True)
