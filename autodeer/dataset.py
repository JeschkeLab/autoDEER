
import numpy as np
import matplotlib.pyplot as plt
import os
import uuid
import json
import base64
import numbers
from scipy.signal import hilbert
from scipy.io import savemat
from autodeer.utils import autoEPRDecoder
from autodeer import __version__
from autodeer.pulses import Pulse,Detection
from autodeer.classes import Parameter, Interface
from autodeer.sequences import Sequence
import autodeer.pulses as ad_pulses
import autodeer.sequences as ad_seqs
import xarray as xr
from deerlab import correctphase
from deerlab import deerload

# =============================================================================
class Dataset:
    """
    Represents an experimental dataset.
    """
    def __init__(self, axes: np.ndarray, data: np.ndarray, params: dict = None,
                 scans: np.ndarray = None, force_complex:bool = True,
                 num_scans = 0, **kwargs) -> None:
        """
        Parameters
        ----------
        axes : list[np.ndarray]
            An array of vectors containing the axes of the dataset
        data : np.ndarray
            The data either as an array of n-dimensions, containing either 
            float or complex values.
        params : dict, optional
            A dictionary of experimental parameters, by default None
        """
        if type(axes) != list:
            self.axes = [axes]
        else:
            self.axes = axes
        self.data = data
        self.params = params
        if type(self.axes) is np.ndarray:
            self.dims = self.axes.ndim
        else:
            self.dims = len(self.axes)
        self.scans = scans
        self.num_scans = Parameter(name='num_scans', value=num_scans)
        if force_complex:
            if not np.iscomplexobj(self.data):
                self.data = hilbert(self.data)
        pass

    def plot(self, label=None, **kwargs) -> plt.figure:
        """
        Produces a standard quick graph of the data. This is a line plot for 
        1D data and a heatmap for 2D data.

        Parameters
        ----------
        label : list[str], optional
            A list contains labels. [x_label,z_label], by default None
        
        Optional Kwargs
        ----------
        lines : list[str], optional
            A list of lines to plot. A selection of ["real", "imag","abs"]
            by default all options are plotted.

        Returns
        -------
        plt.figure
            _description_
        """

        if "lines" in kwargs:
            if 'abs' in kwargs["lines"]:
                abs = True
            if 'imag' in kwargs["lines"]:
                imag = True
            if 'real' in kwargs["lines"]:
                real = True
        else:
            abs=True
            real=True
            imag=True
        if self.dims == 1:
            fig, ax = plt.subplots()
            if abs:
                ax.plot(self.axes[0], np.abs(self.data), label='abs')
            if real:
                ax.plot(self.axes[0], np.real(self.data), label='real')
            if imag:
                ax.plot(self.axes[0], np.imag(self.data), label='imag')
            ax.legend()
            ax.set_ylabel('Signal')
            if label is not None:
                ax.set_xlabel(label[0])
        else:
            raise RuntimeError("Only Single dimension data is supported")

        return fig
        
    def save(self, file: str) -> None:
        """
        Saves the dataset in a variety of commonly used data formats based of
        the file extension. 

        Extension options: [".mat",".np",".txt"]

        Parameters
        ----------
        file : str
            The file name or path to save to. 
        """
        filename, file_ext = os.path.splittext(file)
        
        if file_ext == ".mat":  # Save as old-style .mat file
            save_dict = {
                "dta": self.data,
                "axes": self.axes,
                "params": self.params
                }
            savemat(filename, save_dict) 
        elif file_ext == ".np":  # Save as numpy file
            if self.dims == 1:
                save_array = np.vstack([self.axes, self.data])
                np.save(filename, save_array)

        elif file_ext == ".txt":  # Save as text file
            if self.dims == 1:
                save_array = np.vstack([self.axes, self.data])
                np.savetxt(filename, save_array, delimiter=",")

    def add_variable(self, param):
        setattr(self, param.name, param)


    def _to_dict(self):
        to_return = {"version": __version__, "type": "Sequence"}

        for key, var in vars(self).items():
            if isinstance(var, Parameter):
                to_return[key] = var._to_dict()
            if key == "pulses":
                new_list = []
                for pulse in var:
                    new_list.append(pulse._to_dict())
                to_return[key] = new_list
            else:
                to_return[key] = var

        return to_return

    def _to_json(self):
        class autoEPREncoder(json.JSONEncoder):
            def default(self, obj):
                if isinstance(obj, np.ndarray):
                    if (len(obj) > 0 ) and isinstance(obj[0], str):
                        return list(obj)
                    data = np.ascontiguousarray(obj.data)
                    data_b64 = base64.b64encode(data)
                    return dict(__ndarray__=str(data_b64),
                                dtype=str(obj.dtype),
                                shape=obj.shape)
                if isinstance(obj, complex):
                    return str(obj)
                if isinstance(obj, numbers.Number):
                    return float(obj)
                if isinstance(obj, uuid.UUID):
                    return_dict = {"__uuid__": str(obj)}
                    return return_dict
                if isinstance(obj, Parameter):
                    return obj._to_dict()
                if isinstance(obj, Pulse):
                    return obj._to_dict()
                if isinstance(obj, Sequence):
                    return obj._to_dict()
                else:
                    return json.JSONEncoder.default(self, obj)
        
        return json.dumps(self._to_dict(), cls=autoEPREncoder, indent=4)

    def save(self, filename):
        """Save the sequence to a JSON file.

        Parameters
        ----------
        filename : str
            Path to the JSON file.

        Returns
        -------
        None

        Raises
        ------
        TypeError
            If the object cannot be serialized to JSON.

        Example
        -------
        >>> obj = Sequence()
        >>> obj.save("my_sequence.json")
        """
        
                
        with open(filename, "w") as f:
           f.write(self._to_json())

    @classmethod
    def _from_dict(cls, dct):
        new_dataset = cls(
            axes=dct['axes'], data=dct['data'],
            params=dct['params'],scans = dct['scans'])
        
        for key, var in dct.items():

            if isinstance(var, dict) and ("type" in var):
                
                if var["type"] == "Parameter":
                    setattr(new_dataset, key, Parameter._from_dict(var))
                
                elif var["type"] == "Pulse":
                    try:
                        setattr(new_dataset, key, getattr(ad_pulses, var['subclass'])._from_dict(var))
                    except:
                        setattr(new_dataset, key, Pulse._from_dict(var))

                elif var["type"] == "Sequence":
                    try:
                        setattr(new_dataset, key, getattr(ad_seqs, var['subclass'])._from_dict(var))
                    except:
                        setattr(new_dataset, key, Sequence._from_dict(var))
                
                # elif var["type"] == "Interface":
                    # try:
                    #     setattr(new_dataset, key, getattr(hware, var['subclass'])._from_dict(var))
                #     setattr(new_dataset, key, Interface._from_dict(var))
            
            elif key == "version" or key == "type" or key == "subclass":
                continue
            
            else:
                setattr(new_dataset, key, var)
        
        return new_dataset
    
    @classmethod
    def _from_json(cls, JSONstring):
        dct = json.loads(JSONstring, object_hook=autoEPRDecoder)
        return cls._from_dict(dct)
    
    @classmethod
    def load(cls, filename):
        """Load an object from a JSON file.

        Parameters
        ----------
        filename : str
            Path to the JSON file.

        Returns
        -------
        obj : Dataset
            The Dataset loaded from the JSON file.

        Raises
        ------
        FileNotFoundError
            If the file does not exist.

        Example
        -------
        >>> obj = Dataset.load("my_dataset.json")
        """
        with open(filename, "r") as f:
           file_buffer = f.read()
        return cls._from_json(file_buffer)
    
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
        if not isinstance(param, Parameter):
            continue
        if param.axis == []:
            fixed_param[param_name] = param.value
        

    for i,pulse in enumerate(sequence.pulses):
        if isinstance(pulse, Detection):
            type="det"
        else:
            type="pulse"
        for param_name in pulse.__dict__:
            param = pulse.__dict__[param_name]
            if not isinstance(param, Parameter):
                continue
            if param.axis == []:
                fixed_param[f"{type}{i}_{param_name}"] = param.value
        
    return fixed_param

def create_dataset_from_sequence(data, sequence: Sequence):
    ndims = data.ndim
    default_labels = ['X','Y','Z','T']
    dims = default_labels[:ndims]
    axes = get_all_axes(sequence)
    coords = {a:(default_labels[b['ax_id']],b['axis']) for a,b in axes.items()}
    attr = get_all_fixed_param(sequence)
    
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
    
    return xr.DataArray(data, dims=dims, coords=coords, attrs=attr)

@xr.register_dataarray_accessor("epr")
class EPRAccessor:

    def __init__(self, xarray_obj):
        self._obj = xarray_obj

    @property
    def save(self, filename,type='netCDF'):

        if type == 'netCDF':
            self._obj.to_netcdf(f"{filename}.epr")

    @property
    def correctphase(self, full_output=False, offset=False):


        if np.iscomplexobj(self._obj.data):
            self._obj.data = correctphase(self._obj.data, full_output=full_output, offset=offset)
        else:
            UserWarning("Data is not complex, phase correction not applied")
        return self._obj
    
    @property
    def SNR(self):
        from deerlab import der_snr
        norm_data = self._obj.data / np.abs(self._obj.data).max()
        return 1/der_snr(norm_data)  

        
        
        

