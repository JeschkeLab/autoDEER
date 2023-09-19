from autodeer.classes import  Interface, Parameter
from autodeer.dataset import  Dataset
from autodeer.pulses import Pulse, RectPulse, ChirpPulse, HSPulse, Delay, Detection
from autodeer.sequences import Sequence, HahnEchoSequence, DEERSequence, FieldSweepSequence, ResonatorProfileSequence
from autodeer.FieldSweep import create_Nmodel

import numpy as np
import deerlab as dl
import time

rng = np.random.default_rng(12345)
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


    def __init__(self, speedup=100) -> None:
        self.state = False
        self.speedup = speedup
        self.pulses = {}
        self.start_time = 0

        # Create virtual mode

        def lorenz_fcn(x, centre, sigma):
            y = (0.5*sigma)/((x-centre)**2 + (0.5*sigma)**2)
            return y

        mode = lambda x: lorenz_fcn(x, 34.0, 34.0/60)
        x = np.linspace(33,35)
        scale = 75/mode(x).max()
        self.mode = lambda x: lorenz_fcn(x, 34.0, 34.0/60) * scale
        super().__init__()

    def launch(self, sequence, savename: str, **kwargs):
        self.state = True
        self.sequence = sequence
        self.start_time = time.time()
        return super().launch(sequence, savename)
    
    def acquire_dataset(self,**kwargs) -> Dataset:
        if isinstance(self.sequence, DEERSequence):
            if self.sequence.t.is_static():
                axes, data = _simulate_CP(self.sequence)
            else:
                axes, data =_simulate_deer(self.sequence)
        elif isinstance(self.sequence, FieldSweepSequence):
            axes, data = _simulate_field_sweep(self.sequence)
        elif isinstance(self.sequence,ResonatorProfileSequence):
            axes, data = _similate_respro(self.sequence,self.mode)
        
        data = add_noise(data, 1e-2)
        dset = Dataset(axes, data, scans=1)
        dset.LO = self.sequence.LO
        dset.sequence = self.sequence
        
        
        return dset
    
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

        pulse.scale = Parameter('scale',0.5,unit=None,description='The amplitude of the pulse 0-1')
        return pulse
            
    def isrunning(self) -> bool:
        current_time = time.time()
        if current_time - self.start_time > (self.sequence._estimate_time() / self.speedup):
            self.state = False

        return self.state
    
    def terminate(self) -> None:
        self.state = False
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
    data = func(xaxis,1,0.0000008,1.8)
    data = add_phaseshift(data, 0.05)
    return xaxis, data


def _similate_respro(sequence, mode):

    damped_oscilations = lambda x, f, c: np.sin(2*np.pi*f*x) * np.exp(-c*x)
    damped_oscilations_vec = np.vectorize(damped_oscilations)
    LO_axis = sequence.LO.value + sequence.LO.axis[0]['axis']
    LO_len = LO_axis.shape[0]
    tp_x = val_in_ns(sequence.pulses[0].tp)
    tp_len = tp_x.shape[0]
    nut_freqs = mode(LO_axis)

    damped_oscilations_vec
    data = damped_oscilations_vec(tp_x.reshape(tp_len,1),nut_freqs.reshape(1,LO_len)*1e-3,0.06)
    return [tp_x, LO_axis], data
