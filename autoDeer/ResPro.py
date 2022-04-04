## This is a script for function that run a resonator profile experiment and analyse the data. 

import imp
import time
import numpy as np
import scipy.fft as fft
import importlib
from deerlab import deerload
import re
import os
from scipy.interpolate import pchip_interpolate
from scipy.signal import hilbert,windows
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

MODULE_DIR = importlib.util.find_spec('autoDeer').submodule_search_locations

def run_nutation(api,pulse_lengths:dict={"p0":16,"p1":16,"p2":200},delays:dict={"d1":400,"d7":1000},steps,avgs):
    
    # exp_file = '/home/xuser/Desktop/huka/autoDeer/autoDeer/PulseSpel/autoDEER_4p.exp'
    exp_file = MODULE_DIR + "/PulseSpel/HUKA_DEER_AWG.exp"
    # def_file = '/home/xuser/Desktop/huka/autoDeer/autoDeer/PulseSpel/autoDEER_4p.def'
    def_file = MODULE_DIR+"/PulseSpel/HUKA_DEER_AWG.def"

    # 
    api.set_ReplaceMode(False) #Turn replace mode off
    api.set_set_PhaseCycle(True)
    

    api.set_PulseSpel_exp_filepath(exp_file)
    api.set_PulseSpel_def_filepath(def_file)
    api.compile_PulseSpel_prg()
    api.compile_PulseSpel_def()         
    
    # Set Pulse Lengths
    api.set_PulseSpel_var("p0",pulse_lengths[1])
    api.set_PulseSpel_var("p1",pulse_lengths[0])
    api.set_PulseSpel_var("p2",pulse_lengths[2])
    
    # Set delays
    api.set_PulseSpel_var("d0",delays[0])
    api.set_PulseSpel_var("d1",delays[1])
    api.set_PulseSpel_var("d2",delays[2])
    # Set Steps
    api.set_PulseSpel_var("d30",steps[0])
    api.set_PulseSpel_var("d31",steps[1])
    api.set_PulseSpel_var("d29",steps[2])
    # Set counters
    api.set_PulseSpel_var("h",avgs[0])
    api.set_PulseSpel_var("n",avgs[1])
    api.set_PulseSpel_var("m",avgs[2])
    
    api.set_PulseSpel_experiment("AWG seq inv (p2)")
    api.set_PulseSpel_phase_cycling("AWG +-<x> obs")
    api.set_Acquistion_mode(1)
    
    # Compile Defs and Program
    api.compile_PulseSpel_prg()
    api.compile_PulseSpel_def()  

    # Run Experiment
    api.run_exp()
    time.sleep(1)
    return 1


def get_nutations(api,nu,field,step,nx:int=128):

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

class resonatorProfile:

    def __init__(self,IF_freq=1.25) -> None:
        self.IF_freq =  IF_freq
        pass

    def calc_res_prof(self,nuts:np.ndarray,dt:float=1,throwlist:list=None,freq_lim:tuple=None):

        if nuts.ndim != 2:
            raise ValueError("nuts must be a 2D array where collumn 1 is the frequency")
        n_files = np.shape(nuts)[0]
        nx = np.shape(nuts)[1] - 1

        nut_data = nuts[:,1:]
        nu1_cutoff = 5*1e-3

        fax = fft.fftfreq(nx*2*10,d=dt)
        fax = fft.fftshift(fax)

        fnut = fft.fft(nut_data,nx*2*10,axis=1)
        fnut = fft.fftshift(fnut,axes=1)

        nuoff = np.argwhere(fax > nu1_cutoff)
        nu1_id = np.argmax(np.abs(fnut[:,nuoff]),axis=1)
        nu1 = abs(fax[nu1_id + nuoff[0]])

        data = np.transpose([abs(nuts[:,0]),abs(fax[nu1_id + nuoff[0]]).reshape(n_files)]) 

        throw_mask = np.zeros(0,dtype=np.int16)
        if freq_lim != None:
            throw_mask = np.squeeze(np.argwhere(data[:,0]<freq_lim[0]))
            throw_mask = np.append(throw_mask,np.squeeze(np.argwhere(data[:,0]>freq_lim[1])))

        if throwlist != None:
            for freq in throwlist:
                throw_mask = np.append(throw_mask,np.squeeze(np.argwhere(abs(data[:,0] - freq)<0.001)))

        data = np.delete(data,throw_mask,axis=0)

        self.res_prof = data
    
        return self.res_prof

    def import_from_bruker(self,folder_path,reg_exp="nut_[0-9]+.DSC"):
        
        files = os.listdir(folder_path)
        r = re.compile(reg_exp)
        dsc_files = list(filter(r.match,files))
        n_files = len(dsc_files)
        
        # Load one file to get the time axis and length
        t,V = deerload(folder_path + files[0], False, False)
        t_axis = t
        nx = len(V)

        nut_data = np.zeros((n_files,nx+1),dtype=np.complex64)
        
        for i,files in enumerate(dsc_files):
            t,V,params = deerload(folder_path + files, False, True)
            if not np.array_equal(t,t_axis):
                raise Exception("all files must have the same time axis")
            spec_freq = re.findall("[0-9]+.[0-9]*",params['DSL']['freqCounter']['FrequencyMon'])[0]
            nut_data[i,0] = spec_freq
            nut_data[i,1:] = V

        # sort nut_data
        nut_data_sorted =nut_data[nut_data[:,0].argsort(),:]

        return t_axis, nut_data_sorted

    def calculate_shape(self,data,params=None):
        if params == None:
            params = self.params
        
        A = params[0]
        self.q = params[1]
        self.fc = params[2]

        nu_x = data[:,0]
        nu_y = data[:,1] * 1e3 # Convert from GHz to MHz

        nu_y = nu_y/A

        # Build domains of filters
        srate = 36 * 3 * 4
        # This comes from Q/W = Bandwidth. So we want our frequency range to be 10 * Bandwidth.
        # w = 2 *pi * 1.25(IF freq)
        npts = int(2 ** (np.ceil(np.log2(10 * 50 / (2 * np.pi * 1.25 )* srate ))+ 1))
        frq = np.linspace(-1/2,1/2,npts,endpoint=False) * srate
        habs = np.zeros(frq.shape)

        # Interpolate the resonator profile onto the new frequency axis
        int_idx = np.squeeze(np.argwhere((nu_x[0]< frq) & (frq < nu_x[-1]))) # Where to interpolate
        habs[int_idx] = pchip_interpolate(nu_x,nu_y,frq[int_idx])
        habs_o = habs
        habs = habs/np.max(habs)
        self.nu_max = habs_o[frq>33].max()

        # Lorentzian model
        z = self.q * (frq**2 -self.fc**2)/(frq*self.fc)
        model = (1-z*1j )/ (1 + z**2) # Lorentzian = 1-zj/(1+z^2)
        model[np.isnan(model)] = 0 
        l_abs = np.abs(model)  # Model fit
        l_abs /= np.nanmax(l_abs) # Normalization
        zrs = np.squeeze(np.argwhere(l_abs==0))
        l_abs[zrs] = l_abs[zrs+1]

        l_h = l_abs[len(l_abs)//2+1:-1] # Higher frequencies, above zero
        l_abs_ds = np.concatenate([np.flip(l_h), l_h]) # Make symetrical about zero
        l_pha = np.imag(hilbert(-np.log(l_abs_ds))) # Get phase from the hilbert transform

        # Adding the exp resonantor profile into a lorentzian function
        fdec = 0.1/np.log(2)
        n_dec = int(np.ceil(fdec/(frq[1]-frq[0])))
        n_cyc = 10

        # left side
        habs[0:int_idx[0]] = l_abs[0:int_idx[0]] # put into lorentzian

        # smooth with exp function
        s_range = np.arange(int_idx[0]-n_cyc*n_dec,int_idx[0]+1)
        habs[s_range] = l_abs[s_range] + np.flip((habs[int_idx[0]] - l_abs[int_idx[0]])*np.exp(-1 * np.arange(0,n_cyc+1/n_dec,1/n_dec)))

        # #smooth the edges with a window
        # w_pts = 10
        # win = windows.gaussian(w_pts,2.5) / np.sum(windows.gaussian(w_pts,2.5))
        # w_range = np.arange(int_idx[0]-2*w_pts,int_idx[0] + 2*w_pts)
        # sm = np.convolve(habs[w_range],win,mode='valid')
        # ins_range = np.arange(w_range[0] + np.ceil(w_pts/2),w_range[0] - np.floor(w_pts/2))

        # Repeat for the right hand side
        habs[int_idx[-1]:-1] = l_abs[int_idx[-1]:-1]

        # Smooth with exp function
        s_range = np.arange(int_idx[-1],int_idx[-1]+n_cyc*n_dec+1)
        habs[s_range] = l_abs[s_range] + (habs[int_idx[-1]-1] - l_abs[int_idx[-1]-1])*np.exp(-1 * np.arange(0,n_cyc+1/n_dec,1/n_dec))

        h_pha = np.imag(hilbert(-np.log(habs)))

        self.habs = habs
        self.labs = l_abs
        self.l_pha = l_pha
        self.h_pha = h_pha
        self.frq = frq

    def phase_plots(self,nu_lims:tuple=None,xaxis = 'Bridge'):
        
        if xaxis.lower() == 'bridge':
            frq_axis = self.frq
        elif xaxis.lower() == 'IF':
            frq_axis = self.frq - 32 # Check this number

        if nu_lims != None:
            fmin = nu_lims[0]
            fmax = nu_lims[1]
        else:
            fmin = 0
            fmax = 40 
        
        fig, ax = plt.subplots()
        ax.plot(frq_axis,self.habs,label="Interpolated fit")
        ax.plot(frq_axis,self.labs,label="Lorentzian Model")
        ax.set_xlabel('Frequency GHz')
        ax.set_ylabel('Normalised Nutation Freq.')
        ax.legend()
        ax.set_title('Resontor Profile')
        ax.set_xlim((fmin,fmax))

        return fig

    def res_prof_plot(self,nu_lims:tuple=None,xaxis = 'Bridge'):


        if xaxis.lower() == 'bridge':
            frq_axis = self.frq
        elif xaxis.lower() == 'IF':
            frq_axis = self.frq - 32 # Check this number

        if nu_lims != None:
            fmin = nu_lims[0]
            fmax = nu_lims[1]
        else:
            fmin = 0
            fmax = 40 
        
        fig, ax = plt.subplots()
        ax.plot(frq_axis,self.habs,label="Interpolated fit")
        ax.plot(frq_axis,self.labs,label="Lorentzian Model")
        ax.set_xlabel('Frequency GHz')
        ax.set_ylabel('Normalised Nutation Freq.')
        ax.legend()
        ax.set_title('Resonator Profile')
        ax.set_xlim((fmin,fmax))

        return fig

    
    def autofit(self,nu_lims:tuple=None):
        def lorentz(f,a,q,fc):
            z = q * (f**2 -fc**2)/(f*fc)
            model = a*(1-z*1j )/ (1 + z**2) # Lorentzian = 1-zj/(1+z^2)
            return np.abs(model)

        if nu_lims != None:
            fmin = nu_lims[0]
            fmax = nu_lims[1]
        else:
            fmin = 33
            fmax = 35 

        nu_x = self.res_prof[:,0]
        nu_y = self.res_prof[:,1] *1e3
        nu_y_norm = nu_y/nu_y.max() 

        res, pcov = curve_fit(lorentz,nu_x,nu_y_norm,[1,100,34.5],bounds=([0.7,0,33],[1.3,1000,35]))
        std = np.sqrt(np.diag(pcov))

        self.params = res
        self.q = self.params[1]
        self.fc =self.params[2]
        self.params_std = std


    def calc_IF_prof(self,IF:float,bridgeFreq:float,nu_lims = [33,35],type = 'fit'):
        int_idx = np.squeeze(np.argwhere((nu_lims[0]< self.frq) & (self.frq < nu_lims[-1]))) # Where to interpolate
        self.IF = self.frq[int_idx] - bridgeFreq + IF
        if type == 'fit':
            self.IF_rp = self.labs[int_idx]
        elif type == 'raw':
            self.IF_rp = self.habs[int_idx]


        






