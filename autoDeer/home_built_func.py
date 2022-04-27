import imp
import os,re
import numpy as np
from scipy.io import loadmat
import scipy.fft as fft
import scipy.signal as sig


def find_AWG_deer_files(folder,start_time,experiment="4pdeer",end_time=None):
    files_list = os.listdir(folder)
    include_files_list = []
    re_mask = r"([0-9]{8})_([0-9]{4})_"+experiment+".*.mat"
    for file in files_list:
        if re.match(re_mask,file):
            include = True
            extract = re.findall(re_mask,file)
            if end_time != None:
                if int(extract[0][0]) > end_time[0]:
                    include = False
                elif (int(extract[0][0]) == end_time[0])& (int(extract[0][1]) > end_time[1]):
                    include = False
            elif include:
                if int(extract[0][0]) < start_time[0]:
                    include = False
                elif (int(extract[0][0]) == start_time[0]) & (int(extract[0][1]) < start_time[1]):
                    include = False
            if include:
                include_files_list.append(file)

    return include_files_list


def calc_percieved_freq(f_sample,f_signal):
    return np.abs(f_signal-f_sample*np.round(f_signal/f_sample))

def apply_match(data,freq,filter='rect',sampling_rate=2,integrate=False):
    # Filter is always applied in frequency domain
    filter='rect'
    integrate=True
    data = sig.hilbert(data)
    ft_data = fft.fftshift(fft.fft(data))
    ft_axis = fft.fftshift(fft.fftfreq(ft_data.shape[0],1/sampling_rate))
    t_axis = np.linspace(0,data.shape[0]*1/sampling_rate,data.shape[0])
    filt_pos = np.argmin(np.abs(freq-ft_axis))
    if filter == 'rect':
        def rect(center,width,length):
            data = np.ones(length)
            data[0:center-width] = 0
            data[center+width:-1] = 0
            return data
        filt = rect(filt_pos,10,ft_axis.shape[0])
    ft_filt_data = ft_data * filt
    time_data= fft.ifft(fft.fftshift(ft_data))
    time_filt_data = fft.ifft(fft.fftshift(ft_filt_data))
    time_filt_data_dc = time_filt_data*np.exp(-1j*2*np.pi*freq*t_axis)
    time_data_dc = time_data*np.exp(-1j*2*np.pi*freq*t_axis)

    if integrate:
        int_width = 48*sampling_rate
        centre = 175
        return np.trapz(time_filt_data_dc[centre-int_width//2:centre+int_width//2])
        # return np.trapz(time_data_dc)
    else:
        return time_data_dc


def build_trace(filepath,files_list):
    data_1 = loadmat(filepath+files_list[0],simplify_cells=True)
    t_dim = data_1['deer']['parvars'][1]['dim']
    deadtime = 2*data_1['deer']['events'][1]['t']-data_1['deer']['events'][2]['t']
    tau1 = data_1['deer']['events'][1]['t']
    tau2 = (data_1['deer']['events'][3]['t'] - data_1['deer']['events'][1]['t']) - tau1
    t = data_1['deer']['parvars'][1]['axis'].astype(np.int64)
    t=t - data_1['deer']['events'][2]['t'] - deadtime
    t = t/1000
    trace = np.zeros(t_dim,dtype=np.complex128)

    for file in files_list:
        trace += get_deer(filepath+file)
    
    params = {'tau1':tau1,'tau2':tau2,'deadtime':deadtime}
    return t,trace,params

def get_deer(filepath):
    all_data = loadmat(filepath,simplify_cells=True)
    data = all_data['dta']
    sampling_rate = 2
    nu_obs = all_data['deer']['events'][0]['pulsedef']['nu_init']
    freq = calc_percieved_freq(sampling_rate,nu_obs)
    trace = -1* np.apply_along_axis(apply_match,0,data,freq)
    return trace    