from autodeer.classes import  Interface, Parameter
from autodeer.dataset import  Dataset
from autodeer.pulses import Pulse, RectPulse, ChirpPulse, HSPulse, Delay, Detection
from autodeer.sequences import Sequence, HahnEchoSequence, DEERSequence, FieldSweepSequence, ResonatorProfileSequence
from autodeer.FieldSweep import create_Nmodel

import numpy as np
import deerlab as dl


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

class dummyInterface(Interface):


    def __init__(self, speedup=100) -> None:
        self.state = False
        self.speedup = speedup


        # Create virtual mode

        def lorenz_fcn(x, centre, sigma):
            y = (0.5*sigma)/((x-centre)**2 + (0.5*sigma)**2)
            return y

        mode = lambda x: lorenz_fcn(x, 34.0, 34.0/60)
        x = np.linspace(33,35)
        scale = 75/mode(x).max()
        self.mode = lambda x: lorenz_fcn(x, 34.0, 34.0/60) * scale
        super().__init__()

    def launch(self, sequence, savename: str):
        self.state = True
        self.sequence = sequence
        return super().launch(sequence, savename)
    
    def acquire_dataset(self) -> Dataset:
        if isinstance(self.sequence, DEERSequence):
            if self.sequence.t.is_static():
                axes, data = _simulate_CP(self.sequence)
            else:
                axes, data =_simulate_deer(self.sequence)
        elif isinstance(self.sequence, FieldSweepSequence):
            axes, data = _simulate_field_sweep(self.sequence)
        elif isinstance(self.sequence,ResonatorProfileSequence):
            axes, data = _similate_respro(self.sequence,self.mode)


        dset = Dataset(axes, data, scans=1)
        dset.LO = self.sequence.LO
        dset.sequence = self.sequence
        
        
        return dset
            
    def isrunning(self) -> bool:
        return self.state
    
    def terminate(self) -> None:
        self.state = False
        return super().terminate()
    

def _simulate_field_sweep(sequence):
    # Assuming a Nitroxide sample
    Vmodel = create_Nmodel(sequence.LO.value *1e3)
    axis = sequence.B.value + sequence.B.axis[0]['axis']
    axis = axis * 0.1
    Boffset=0
    gy = 2.0061
    gz = 2.0021
    axy = 0.488
    az = 3.66
    GB = 0.45
    scale=1

    data = Vmodel(axis,Boffset,gy,gz,axy,az,GB,scale)
    return axis,data.astype(np.complex128)


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
        mod_depths = {"lam1":0.4, "lam2":0.05, "lam3":0.2, "lam4":0.05, "lam5":0.1}

    elif exp_type == "3pDEER":
        experimentInfo = dl.ex_3pdeer(tau=tau1,pathways=[1,2])
        reftimes = dict(zip(["reftime1","reftime2"],experimentInfo.reftimes(tau1,)))
        mod_depths = {"lam1":0.6, "lam2":0.1}


    r = np.linspace(0.5,10,100)
    rmean = 3.5
    rstd = 1.5
    conc = 50

    Vmodel = dl.dipolarmodel(t,r,Pmodel=dl.dd_gauss, experiment=experimentInfo)
    Vsim = Vmodel(mean=rmean, std=rstd, conc=conc, scale=1, **reftimes, **mod_depths)
    return t, Vsim

def _simulate_CP(sequence):


    func = lambda x, a, b, e: a*np.exp(-b*x**e)
    xaxis = val_in_ns(sequence.tau2)
    data = func(xaxis,1,0.0000008,1.8)
    return xaxis, 


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
