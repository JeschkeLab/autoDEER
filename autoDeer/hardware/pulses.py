import numpy as np
from autoDeer.ResPro import resonatorProfile
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import pchip_interpolate

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
