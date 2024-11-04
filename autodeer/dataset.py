
import numpy as np
import matplotlib.pyplot as plt
import re
import os
import uuid
import json
import base64
import numbers
from scipy.signal import hilbert
from scipy.io import savemat
import scipy.fft as fft
from autodeer.utils import autoEPRDecoder
from autodeer import __version__
from autodeer.pulses import Pulse,Detection
from autodeer.classes import Parameter, Interface
from autodeer.sequences import Sequence, DEERSequence
import autodeer.pulses as ad_pulses
import autodeer.sequences as ad_seqs
import xarray as xr
from deerlab import correctphase
from deerlab import deerload
import copy

# =============================================================================

def get_all_axes(sequence):
    # loop over all parameters and get the axis
    axes = {}
    for param_name in sequence.__dict__:
        param = sequence.__dict__[param_name]
        if not isinstance(param, Parameter):
            continue
        if param.axis == []:
            continue
        
        for axis in param.axis:
            if axis['uuid'] in sequence.axes_uuid and not axis['uuid'] in sequence.reduce_uuid:
                ax_id = sequence.axes_uuid.index(axis['uuid'])
                if param.unit == 'ns':
                    convert = 1e-3
                elif param.unit == 'us':
                    convert = 1
                elif param.unit == 'ms':
                    convert = 1e3
                else:
                    convert = 1
                axes[param_name] = {'axis': axis['axis']*convert + param.value*convert, 'ax_id':ax_id}
    
    for i,pulse in enumerate(sequence.pulses):
        for param_name in pulse.__dict__:
            param = pulse.__dict__[param_name]
            if not isinstance(param, Parameter):
                continue
            if param.axis == []:
                continue
            for axis in param.axis:
                if axis['uuid'] in sequence.axes_uuid and not axis['uuid'] in sequence.reduce_uuid:
                    ax_id = sequence.axes_uuid.index(axis['uuid'])
                    axes[f"pulse{i}_{param_name}"] = {'axis': axis['axis'] + param.value, 'ax_id':ax_id}

    return axes

def get_all_fixed_param(sequence):

    fixed_param = {}
    for param_name in sequence.__dict__:
        param = sequence.__dict__[param_name]
        if param_name == 'name':
            if param is not None:
                fixed_param[f"seq_{param_name}"] = param
                continue 
        elif not isinstance(param, Parameter):
            continue
        else:
            if (param.axis == []) and (param.value is not None):
                fixed_param[param_name] = param.value
            elif param.axis[0]['uuid'] in sequence.reduce_uuid:
                fixed_param[param_name] = param.value
        

    for i,pulse in enumerate(sequence.pulses):
        if isinstance(pulse, Detection):
            type="det"
        else:
            type="pulse"
        for param_name in pulse.__dict__:
            param = pulse.__dict__[param_name]
            if param_name == 'name':
                if param is not None:
                    fixed_param[f"{type}{i}_{param_name}"] = param
                    continue
            elif not isinstance(param, Parameter):
                continue
            else:
                if (param.axis == []) and (param.value is not None):
                    fixed_param[f"{type}{i}_{param_name}"] = param.value

    fixed_param['nPcyc'] = sequence.pcyc_dets.shape[0]
    if isinstance(sequence,DEERSequence):
        fixed_param['pcyc_name'] = sequence.pcyc_name
    return fixed_param

def create_dataset_from_sequence(data, sequence: Sequence,extra_params={}):
    ndims = data.ndim
    default_labels = ['X','Y','Z','T']
    dims = default_labels[:ndims]
    axes = get_all_axes(sequence)
    coords = {a:(default_labels[b['ax_id']],b['axis']) for a,b in axes.items()}
    attr = get_all_fixed_param(sequence)
    attr.update(extra_params)
    
    return xr.DataArray(data, dims=dims, coords=coords,attrs=attr)

def create_dataset_from_axes(data, axes, params: dict = None,axes_labels=None):
    ndims = data.ndim
    if axes_labels is None:
        default_labels = ['X','Y','Z','T']
    elif len(axes_labels) < ndims:
        default_labels = axes_labels + ['X','Y','Z','T']
    else:
        default_labels = axes_labels
    dims = default_labels[:ndims]
    if not isinstance(axes, list):
        axes = [axes]
    coords = {default_labels.pop(0):a for a in axes}
    
    return xr.DataArray(data, dims=dims, coords=coords, attrs=params)

def create_dataset_from_bruker(filepath):
    axes, data, params = deerload(filepath, plot=False, full_output=True)
    if not isinstance(axes, list):
        axes = [axes]
    ndims = data.ndim
    default_labels = ['X','Y','Z','T']
    dims = default_labels[:ndims]
    labels = []
    for i in range(ndims):
        ax_label = default_labels[i]
        axis_string = params['DESC'][f'{ax_label}UNI']
        if "'" in axis_string:
            axis_string = axis_string.replace("'", "")
        if axis_string == 'G':
            labels.append('B')
            axes[i] = axes[i] * 1e3
        elif axis_string == 'ns':
            labels.append('t')
        else:
            labels.append(None)
    # Count occurens of each label
    label_count = {i:labels.count(i) for i in labels}
    # Add a number to each label if count > 1
    for i in range(ndims):
        if label_count[labels[i]] > 1:
            labels[i] = default_labels[i] + labels[i]
    
    coords = {labels[i]:(default_labels[i],a) for i,a in enumerate(axes)}
    attr = {}
    attr['LO'] = float(params['SPL']['MWFQ'])  / 1e9
    attr['B'] = float(params['SPL']['B0VL']) * 1e4
    attr['reptime'] = float(params['DSL']['ftEpr']['ShotRepTime'].replace('us',''))
    attr['nAvgs'] = int(params['DSL']['recorder']['NbScansAcc'])
    attr['shots'] = int(params['DSL']['ftEpr']['ShotsPLoop'])
    
    return xr.DataArray(data, dims=dims, coords=coords, attrs=attr)

@xr.register_dataarray_accessor("epr")
class EPRAccessor:

    def __init__(self, xarray_obj):
        self._obj = xarray_obj

    def save(self, filename,type='netCDF'):

        if type == 'netCDF':
            self._obj.to_netcdf(f"{filename}.h5",engine='h5netcdf',invalid_netcdf=True)

    @property
    def correctphase(self):
        new_obj = copy.deepcopy(self._obj)
        if np.iscomplexobj(self._obj.data):
                corr_data = correctphase(self._obj.data)
                new_obj.data = corr_data
        else:
            UserWarning("Data is not complex, phase correction not applied")
        return new_obj
    
    @property
    def normalise(self):
        self._obj.data = self._obj.data / np.abs(self._obj.data).max()
        return self._obj
    
    @property
    def correctphasefull(self):
        new_obj = copy.deepcopy(self._obj)
        if np.iscomplexobj(self._obj.data):
                Re,Im,_ = correctphase(self._obj.data,full_output=True)
                new_obj.data = Re + 1j*Im
        else:
            UserWarning("Data is not complex, phase correction not applied")
        return new_obj

    @property
    def SNR(self):
        from deerlab import der_snr
        norm_data = self._obj.data / np.abs(self._obj.data).max()
        return 1/der_snr(norm_data)

    @property
    def fft(self):
        new_obj = copy.deepcopy(self._obj)
        new_obj.data = fft.fftshift(fft.fft(self._obj.data))
        new_coords = {}
        for key, coord in new_obj.coords.items():
            new_coords[key] = (coord.dims[0],fft.fftshift(fft.fftfreq(coord.size, coord[1].data-coord[0].data)))
        
        new_obj = new_obj.assign_coords(**new_coords)
        
        return new_obj
    
    
    @property
    def sequence(self):
        dataset_attrs = self._obj.attrs
        dataset_coords = self._obj.coords

        pulses = len([key for key in dataset_attrs.keys() if re.match(r"pulse\d+_name$", key)])

        seq_param_types = ['seq_name','B','LO','reptime','shots','averages','det_window']
        seq_params = {}

        for param_type in seq_param_types:
            if param_type in dataset_attrs:
                seq_params[param_type] = dataset_attrs[param_type]
            elif param_type in dataset_coords:
                coord = dataset_coords[param_type]
                min = coord.min()
                dim = coord.shape[0]
                step = coord[1] - coord[0]

                seq_params[param_type] = Parameter(name = param_type, value = min, dim=dim, step=step)

        seq_params['name'] = seq_params.pop('seq_name')

        pulses_obj = []
        for i in range(pulses):
            pulse_type = dataset_attrs[f"pulse{i}_name"]
            param_types = ['t','tp','freq','flipangle','scale','order1','order2','init_freq','BW','final_freq','beta']
            params = {}

            for param_type in param_types:
                if f"pulse{i}_{param_type}" in dataset_attrs:
                    params[param_type] = dataset_attrs[f"pulse{i}_{param_type}"]
                elif f"pulse{i}_{param_type}" in dataset_coords:
                    coord = dataset_coords[f"pulse{i}_{param_type}"]
                    min = coord.min()
                    dim = coord.shape[0]
                    step = coord[1] - coord[0]

                    params[param_type] = Parameter(name = param_type, value = min, dim=dim, step=step)
            
            try:
                pulse_build = getattr(ad_pulses,pulse_type)
                pulses_obj.append(pulse_build(**params))
            except:
                continue

        sequence = Sequence(**seq_params)
        sequence.pulses = pulses_obj

        return sequence





        
        
        

