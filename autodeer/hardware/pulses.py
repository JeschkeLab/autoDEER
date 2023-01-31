import numpy as np
from autodeer.ResPro import resonatorProfile
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import pchip_interpolate
from scipy.fft import fft,fftshift,fftfreq
import matplotlib.pyplot as plt

def add_modulation(data,freq,phase,sampling_period):
    num_points = len(data)
    t = np.arange(0,num_points) * sampling_period
    return data * np.sin(freq*t+phase)

def gaus_shape(x,amp=1,offset=0,c=10):
    y = amp * np.exp(-(x-offset)**2/c**2)
    return y

def HSorder1(sampling_rate:float,length:int,nu:tuple,resonator:resonatorProfile=None,HSbeta:int=10,HSorder:tuple=[1,6],phase=0,scale=1):


    # Make time axis
    dt = 1/sampling_rate
    n = round(length/dt)
    tp = n*dt
    t = np.linspace(0,tp-dt,n)

    deltaf = nu[1]-nu[0]
    tcent = tp / 2

    # Incoporate AM
    beta = HSbeta
    order = HSorder[0]
    beta_exp = np.log(beta*0.5**(1-order)) / np.log(beta)
    cut = round(len(t)/2)
    AM = np.zeros(n)
    AM[0:cut] = 1/np.cosh(beta**beta_exp * ((t[0:cut] - tcent)/tp)**order)
    order = HSorder[1]
    beta_exp = np.log(beta*0.5**(1-order)) / np.log(beta)
    AM[cut:] = 1/np.cosh(beta**beta_exp * ((t[cut:] - tcent)/tp)**order)

    FM = deltaf * cumulative_trapezoid(AM**2,t,initial=0) / np.trapz(AM**2,t) + nu[0]
    frange = FM
    smoothing = AM

    # Default the resonator profile 
    if resonator == None:
        res_range = [0, 20]
        res_nu1 = [1, 1]
    else:
        res_range = [resonator.IF.min(),resonator.IF.max()]
        res_range = resonator.IF
        res_nu1 = resonator.IF_rp

    sprofile = pchip_interpolate(res_range,np.abs(res_nu1),frange)
    sprofile_orig = sprofile

    sprofile = sprofile * smoothing

    mean_alpha = np.trapz(1 / (sprofile ** 2),frange)/t[-1]
    dt_df = (1 / (sprofile ** 2) )/ mean_alpha
    qref = 2*np.pi / mean_alpha

    # check for nu1 in MHz
    if max(np.abs(res_nu1)) > 1:
        qref = qref * 0.001**2
    
    # Get time tau vs freq from dt/df
    tau = cumulative_trapezoid(dt_df,frange,initial=0)

    # Reparamterise tau(nu) to nu(tau)
    nu = pchip_interpolate(tau,frange,t)

    if (frange[1] - frange[0])<0:
        AM = np.flip(pchip_interpolate(np.flip(frange),np.flip(smoothing),np.flip(nu)))
    else:
        AM = pchip_interpolate(frange,smoothing,nu)

    FM = nu

    real_wave = scale* AM * np.sin(2 * np.pi * cumulative_trapezoid(FM,initial=0) * dt + phase)
    imag_wave = scale* AM * np.cos(2 * np.pi * cumulative_trapezoid(FM,initial=0) * dt + phase)

    return [real_wave,imag_wave]


def rectPulse(sampling_rate:float,length:int,nu,scale=1,phase=0):

    dt = 1/sampling_rate
    n = round(length/dt)    
    tp = n*dt
    t = np.linspace(0,tp-dt,n)

    real_wave = scale * np.sin(2 * np.pi * nu * t + phase)
    imag_wave = scale * np.cos(2 * np.pi * nu * t + phase)

    return [real_wave,imag_wave]


class Pulse:

    def __init__(self,sampling_freq) -> None:
        self.sampling_freq =sampling_freq
        pass

    def calc_exc_profile(self):
        if not hasattr(self,"IQ"):
            raise ValueError("IQ is not defined")
        if not hasattr(self,"sampling_freq"):
            raise ValueError("Sampling Freq is not defined")

        num_points = len(self.IQ)

        self.exc_profile = fftshift(fft(self.IQ))
        self.exc_profile /= self.exc_profile.max()
        self.exc_freqs = fftshift(fftfreq(num_points,1/self.sampling_freq))

        return self.exc_freqs,self.exc_profile

    def plot_exc_profile(self,fieldsweep=None):

        if not hasattr(self,"exc_prof"):
            self.calc_exc_profile()

        fig, ax = plt.subplots()
        ax.plot(self.exc_freqs,np.abs(self.exc_profile),label="Pulse Profile")

        if fieldsweep != None:
            ax.plot(fieldsweep.fs_x,np.abs(fieldsweep.data)/np.abs(fieldsweep.data).max(),label="Field Sweep")
        
        ax.legend()

        return fig

    def pad(self,width,value=0):

        self.IQ = np.pad(self.IQ,width,mode="constant",constant_values=value)
        return self.IQ


class RectPulse(Pulse):

    def __init__(self,sampling_freq,length,nu,scale=1,phase=0) -> None:
        super().__init__(sampling_freq)
        self.sampling_freq = sampling_freq
        self.length = length
        self.nu = nu
        self.scale = scale
        self.phase = phase
        self._build_pulse()


    def _build_pulse(self):

        dt = 1/self.sampling_rate
        n = round(self.length/dt) # Number of points   
        tp = n*dt #This is the actual pulse length
        t = np.linspace(0,tp-dt,n)

        real_wave = self.scale * np.sin(2 * np.pi * self.nu * t + self.phase)
        imag_wave = self.scale * np.cos(2 * np.pi * self.nu * t + self.phase)

        self.IQ = real_wave + 1j*  imag_wave
        return self.IQ


class HSPulse(Pulse):

    def __init__(self,sampling_freq,length:int,nu:tuple,resonator:resonatorProfile=None,HSbeta:int=10,
                HSorder:tuple=[1,6],phase=0,scale=1) -> None:
        super().__init__(sampling_freq)
        self.sampling_freq = sampling_freq
        self.length = length
        self.nu = nu
        self.scale = scale
        self.phase = phase
        self.res_prof = resonator
        self.HSbeta = HSbeta
        self.HSorder = HSorder
        self._build_pulse()

    def _build_pulse(self):
            # Make time axis
        dt = 1/self.sampling_freq
        n = round(self.length/dt)
        tp = n*dt
        t = np.linspace(0,tp-dt,n)

        deltaf = self.nu[1]-self.nu[0]
        tcent = tp / 2

        # Incoporate AM
        beta = self.HSbeta
        order = self.HSorder[0]
        beta_exp = np.log(beta*0.5**(1-order)) / np.log(beta)
        cut = round(len(t)/2)
        AM = np.zeros(n)
        AM[0:cut] = 1/np.cosh(beta**beta_exp * ((t[0:cut] - tcent)/tp)**order)
        order = self.HSorder[1]
        beta_exp = np.log(beta*0.5**(1-order)) / np.log(beta)
        AM[cut:] = 1/np.cosh(beta**beta_exp * ((t[cut:] - tcent)/tp)**order)

        FM = deltaf * cumulative_trapezoid(AM**2,t,initial=0) / np.trapz(AM**2,t) + self.nu[0]
        frange = FM
        smoothing = AM

        # Default the resonator profile 
        if self.res_prof == None:
            res_range = [0, 20]
            res_nu1 = [1, 1]
        else:
            res_range = [self.res_prof.IF.min(),self.res_prof.IF.max()]
            res_range = self.res_prof.IF
            res_nu1 = self.res_prof.IF_rp

        sprofile = pchip_interpolate(res_range,np.abs(res_nu1),frange)
        sprofile_orig = sprofile

        sprofile = sprofile * smoothing

        mean_alpha = np.trapz(1 / (sprofile ** 2),frange)/t[-1]
        dt_df = (1 / (sprofile ** 2) )/ mean_alpha
        qref = 2*np.pi / mean_alpha

        # check for nu1 in MHz
        if max(np.abs(res_nu1)) > 1:
            qref = qref * 0.001**2
        
        # Get time tau vs freq from dt/df
        tau = cumulative_trapezoid(dt_df,frange,initial=0)

        # Reparamterise tau(nu) to nu(tau)
        nu = pchip_interpolate(tau,frange,t)

        if (frange[1] - frange[0])<0:
            AM = np.flip(pchip_interpolate(np.flip(frange),np.flip(smoothing),np.flip(nu)))
        else:
            AM = pchip_interpolate(frange,smoothing,nu)

        FM = nu

        real_wave = self.scale* AM * np.sin(2 * np.pi * cumulative_trapezoid(FM,initial=0) * dt + self.phase)
        imag_wave = self.scale* AM * np.cos(2 * np.pi * cumulative_trapezoid(FM,initial=0) * dt + self.phase)

        self.IQ = real_wave + 1j*imag_wave

