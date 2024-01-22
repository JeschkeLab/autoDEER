from autodeer.classes import  Interface, Parameter
from autodeer.dataset import  create_dataset_from_sequence
from autodeer.pulses import Pulse, RectPulse, ChirpPulse, HSPulse, Delay, Detection
from autodeer.sequences import *
from autodeer.FieldSweep import create_Nmodel
import yaml

import numpy as np
import deerlab as dl
import time
import logging


rng = np.random.default_rng(12345)

hw_log = logging.getLogger('interface.Dummy')

def val_in_us(Param):
        if len(Param.axis) == 0:
            if Param.unit == "us":
                return Param.value
            elif Param.unit == "ns":
                return Param.value / 1e3
        elif len(Param.axis) == 1:
            if Param.unit == "us":
                return Param.tau1.value + Param.axis[0]['axis']
            elif Param.unit == "ns":
                return (Param.value + Param.axis[0]['axis']) / 1e3 

def val_in_ns(Param):
        if len(Param.axis) == 0:
            if Param.unit == "us":
                return Param.value * 1e3
            elif Param.unit == "ns":
                return Param.value 
        elif len(Param.axis) == 1:
            if Param.unit == "us":
                return (Param.tau1.value + Param.axis[0]['axis']) * 1e3
            elif Param.unit == "ns":
                return (Param.value + Param.axis[0]['axis']) 

def add_noise(data, noise_level):
    # Add noise to the data with a given noise level for data that could be either real or complex
    if np.isrealobj(data):
        noise = np.squeeze(rng.normal(0, noise_level, size=(*data.shape,1)).view(np.float64))

    else:
        noise = np.squeeze(rng.normal(0, noise_level, size=(*data.shape,2)).view(np.complex128))
    data = data + noise
    return data

def add_phaseshift(data, phase):
    data = data.astype(np.complex128) * np.exp(-1j*phase*np.pi)
    return data
    

class dummyInterface(Interface):


    def __init__(self,config_file) -> None:
        with open(config_file, mode='r') as file:
            config = yaml.safe_load(file)
            self.config = config
        
        Dummy = config['Spectrometer']['Dummy']
        Bridge = config['Spectrometer']['Bridge']
        resonator_list = list(config['Resonators'].keys())
        self.state = False
        self.speedup = Dummy['speedup']
        self.pulses = {}
        self.start_time = 0
        self.SNR = Dummy['SNR']
        if 'ESEEM_depth' in Dummy.keys():
            self.ESEEM = Dummy['ESEEM_depth']
        else:
            self.ESEEM = 0

        # Create virtual mode
        key1 = resonator_list[0]
        fc = self.config['Resonators'][key1]['Center Freq']
        Q = self.config['Resonators'][key1]['Q']
        def lorenz_fcn(x, centre, sigma):
            y = (0.5*sigma)/((x-centre)**2 + (0.5*sigma)**2)
            return y

        mode = lambda x: lorenz_fcn(x, fc, fc/Q)
        x = np.linspace(33,35)
        scale = 75/mode(x).max()
        self.mode = lambda x: lorenz_fcn(x, fc, fc/Q) * scale
        super().__init__()

    def launch(self, sequence, savename: str, **kwargs):
        hw_log.info(f"Launching {sequence.name} sequence")
        self.state = True
        self.cur_exp = sequence
        self.start_time = time.time()
        return super().launch(sequence, savename)
    
    def acquire_dataset(self,**kwargs):
        hw_log.debug("Acquiring dataset")
        if isinstance(self.cur_exp, DEERSequence):
            if self.cur_exp.t.is_static():
                axes, data = _simulate_CP(self.cur_exp)
            else:
                axes, data =_simulate_deer(self.cur_exp)
        elif isinstance(self.cur_exp, FieldSweepSequence):
            axes, data = _simulate_field_sweep(self.cur_exp)
        elif isinstance(self.cur_exp,ResonatorProfileSequence):
            axes, data = _similate_respro(self.cur_exp,self.mode)
        elif isinstance(self.cur_exp, ReptimeScan):
            axes, data = _simulate_reptimescan(self.cur_exp)
        elif isinstance(self.cur_exp,T2RelaxationSequence):
            axes, data = _simulate_T2(self.cur_exp,self.ESEEM)
        else:
            raise NotImplementedError("Sequence not implemented")

        time_estimate = self.cur_exp._estimate_time()
        if self.speedup != np.inf:
            time_estimate /= self.speedup
            progress = (time.time() - self.start_time) / time_estimate
            if progress > 1:
                progress = 1
            data = add_noise(data, 1/(self.SNR*progress))
        else:
            progress = 1
        scan_num = self.cur_exp.averages.value
        dset = create_dataset_from_sequence(data,self.cur_exp)
        dset.attrs['nAvgs'] = int(scan_num*progress)
        
    
        return super().acquire_dataset(dset)
    
    def tune_rectpulse(self,*,tp, LO, B, reptime,**kwargs):

        rabi_freq = self.mode(LO)
        def Hz2length(x):
            return 1 / ((x/1000)*2)
        rabi_time = Hz2length(rabi_freq)
        if rabi_time > tp:
            p90 = tp
            p180 = tp*2
        else:
            p90 = rabi_time/tp
            p180 = p90*2

        self.pulses[f"p90_{tp}"] = RectPulse(tp=tp, freq=0, flipangle=np.pi/2, scale=p90)
        self.pulses[f"p180_{tp}"] = RectPulse(tp=tp, freq=0, flipangle=np.pi, scale=p180)

        return self.pulses[f"p90_{tp}"], self.pulses[f"p180_{tp}"]
    
    def tune_pulse(self, pulse, mode, LO, B , reptime, shots=400):
        hw_log.debug(f"Tuning {pulse.name} pulse")
        pulse.scale = Parameter('scale',0.5,unit=None,description='The amplitude of the pulse 0-1')
        hw_log.debug(f"Setting {pulse.name} pulse to {pulse.scale.value}")
        return pulse
            
    def isrunning(self) -> bool:
        current_time = time.time()
        if current_time - self.start_time > (self.cur_exp._estimate_time() / self.speedup):
            self.state = False
    
        return self.state
    
    def terminate(self) -> None:
        self.state = False
        hw_log.info("Terminating sequence")
        return super().terminate()
    

def _simulate_field_sweep(sequence):
    # Assuming a Nitroxide sample
    Vmodel = create_Nmodel(sequence.LO.value *1e3)
    axis = sequence.B.value + sequence.B.axis[0]['axis']
    sim_axis = axis * 0.1
    Boffset=0
    gy = 2.0061
    gz = 2.0021
    axy = 0.488
    az = 3.66
    GB = 0.45
    scale=1

    data = Vmodel(sim_axis,Boffset,gy,gz,axy,az,GB,scale)
    data = add_phaseshift(data, 0.05)
    return axis,data


def _simulate_deer(sequence,exp_type=None):

    if sequence.name == "4pDEER":
        exp_type = "4pDEER"
        tau1 = val_in_us(sequence.tau1)
        tau2 = val_in_us(sequence.tau2)
        t = val_in_us(sequence.t)
    elif sequence.name == "5pDEER":
        exp_type = "5pDEER"
        tau1 = val_in_us(sequence.tau1)
        tau2 = val_in_us(sequence.tau2)
        tau3 = val_in_us(sequence.tau3)
        t = val_in_us(sequence.t)
    elif sequence.name == "3pDEER":
        exp_type = "3pDEER"
        tau1 = val_in_us(sequence.tau1)
        t = val_in_us(sequence.t)
    elif sequence.name == "nDEER-CP":
        exp_type = "4pDEER"
        tau1 = val_in_us(sequence.tau1)
        tau2 = val_in_us(sequence.tau2)
        t = val_in_us(sequence.t)

    if exp_type == "4pDEER":
        experimentInfo = dl.ex_4pdeer(tau1=tau1,tau2=tau2,pathways=[1,2,3])
        reftimes = dict(zip(["reftime1","reftime2","reftime3"],experimentInfo.reftimes(tau1,tau2)))
        mod_depths = {"lam1":0.4, "lam2":0.1, "lam3":0.2}
    elif exp_type == "5pDEER":
        experimentInfo = dl.ex_fwd5pdeer(tau1=tau1,tau2=tau2,tau3=tau3,pathways=[1,2,3,4,5])
        reftimes = dict(zip(["reftime1","reftime2","reftime3","reftime4","reftime5"],experimentInfo.reftimes(tau1,tau2,tau3)))
        mod_depths = {"lam1":0.4, "lam2":0.00, "lam3":0.0, "lam4":0.00, "lam5":0.1}

    elif exp_type == "3pDEER":
        experimentInfo = dl.ex_3pdeer(tau=tau1,pathways=[1,2])
        reftimes = dict(zip(["reftime1","reftime2"],experimentInfo.reftimes(tau1,)))
        mod_depths = {"lam1":0.6, "lam2":0.1}


    r = np.linspace(0.5,10,100)
    rmean = 4.5
    rstd = 1.0
    conc = 50

    Vmodel = dl.dipolarmodel(t,r,Pmodel=dl.dd_gauss, experiment=experimentInfo)
    Vsim = Vmodel(mean=rmean, std=rstd, conc=conc, scale=1, **reftimes, **mod_depths)
    # Add phase shift
    Vsim = add_phaseshift(Vsim, 0.05)
    return t, Vsim

def _simulate_CP(sequence):


    func = lambda x, a, b, e: a*np.exp(-b*x**e)
    xaxis = val_in_ns(sequence.tau2)
    data = func(xaxis,1,2e-6,1.8)
    data = add_phaseshift(data, 0.05)
    return xaxis, data

def _simulate_T2(sequence,ESEEM_depth):
    func = lambda x, a, b, e: a*np.exp(-b*x**e)
    xaxis = val_in_ns(sequence.tau)
    data = func(xaxis,1,10e-6,1.6)
    data = add_phaseshift(data, 0.05)
    if ESEEM_depth != 0:
        data *= _gen_ESEEM(xaxis, 7.842, ESEEM_depth)
    return xaxis, data

def _similate_respro(sequence, mode):

    damped_oscilations = lambda x, f, c: np.cos(2*np.pi*f*x) * np.exp(-c*x)
    damped_oscilations_vec = np.vectorize(damped_oscilations)
    LO_axis = sequence.LO.value + sequence.LO.axis[0]['axis']
    LO_len = LO_axis.shape[0]
    tp_x = val_in_ns(sequence.pulses[0].tp)
    tp_len = tp_x.shape[0]
    nut_freqs = mode(LO_axis)

    damped_oscilations_vec
    data = damped_oscilations_vec(tp_x.reshape(tp_len,1),nut_freqs.reshape(1,LO_len)*1e-3,0.06)
    return [tp_x, LO_axis], data

def _simulate_reptimescan(sequence):
    def func(x,T1):
        return 1-np.exp(-x/T1)
    t = sequence.reptime.value + sequence.reptime.axis[0]['axis']
    T1 = 2000 #us

    data = func(t,T1)
    data = add_phaseshift(data, 0.05)
    return t, data

def _simulate_2D_relax(sequence):

    def func(x,y,T2):
        

def _gen_ESEEM(t,freq,depth):
    # Generate an ESEEM modulation
    modulation = np.ones_like(t,dtype=np.float64)
    modulation -= depth *(0.5 + 0.5*np.cos(2*np.pi*t*freq)) + depth * (0.5+0.5*np.cos(2*np.pi*t*freq/2))
    return modulation