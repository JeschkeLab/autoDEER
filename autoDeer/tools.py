import os
import deerlab as dl
import numpy as np
import logging
import h5py
from autodeer.hardware.openepr import dataset, Parameter
from autodeer.home_built_func import uwb_load
from scipy.io import loadmat

log = logging.getLogger('core.Tools')


def eprload(
        path: str, experiment: str = None, type: str = None,
        **kwargs) -> dataset:
    """ A general versions of eprload

    Parameters
    ----------
    path : str
        The file path of the data that should be loaded.
    experiment : str, optional
        _description_, by default None
    type : str, optional
        _description_, by default None

    Returns
    -------
    dataset
        _description_

    Raises
    ------
    ValueError
        _description_
    RuntimeError
        _description_
    """

    if type is None:  # Use the file ending to guess file type
        _, file_extension = os.path.splitext(path)

        if (file_extension == ".DSC") | (file_extension == ".DTA"):
            log.debug('File detected as Bruker')
            type = 'BRUKER'

        elif (file_extension == ".h5") | (file_extension == ".hdf5"):
            log.debug('File detected as HDF5')
            type = 'HDF5'

        elif (file_extension == ".csv") | (file_extension == ".txt"):
            log.debug('File detected as csv or txt file')
            type = 'TXT'

        elif file_extension == '.mat':
            log.debug('File detecetd as Matlab')
            type = 'MAT'
        
        else:
            log.error("Can't detect file type")
            raise ValueError(
                "Can't detect file type. Please choose the correct file or"
                " set type manually \n Valid file types: '.DSC','.DTA','.h5',"
                "'.hdf5','.csv','.txt','.mat'")
    
    if type == 'BRUKER':
        
        if 'full_output' in kwargs:
            full_output = kwargs['full_output']
        else:
            full_output = False

        if full_output is False:    
            t, V = dl.deerload(path, plot=False, full_output=full_output)
            return dataset(t, V)

        else:
            t, V, Params = dl.deerload(
                path, plot=False, full_output=full_output)
            return dataset(t, V, Params)

    elif type == 'HDF5':

        f = h5py.File(path, 'r')
        groups = list(f.keys())
        
        if experiment is None:
            if 'Spectrometer' in groups:
                groups.remove('Spectrometer')
            exp = f[groups[0]]
        else:
            exp = f[experiment]
        
        # attrs = (exp.attrs.keys())
        elements = list(exp.keys())

        if all(item in elements for item in ['Axes', 'Data']):
            t = exp['Axes']
            V = exp['Data']
        else:
            raise RuntimeError(
                "HDF5 File is missing the datasets ('Axes','Data') in group"
                f"{experiment}")

        if 'full_output' in kwargs:
            full_output = kwargs['full_output']
        
            if full_output:
                Params = dict(exp.attrs.items())
                return dataset(t, V, Params)   
            else:
                return dataset(t, V)
        else:
            return dataset(t, V)

    elif type == 'TXT':
        if 'full_output' in kwargs:
            full_output = kwargs['full_output']
            del kwargs['full_output']
            if full_output:
                print("WARNING: Can't get metadata from text file")
        data = np.loadtxt(path, *kwargs)
        return data

    elif type == 'MAT':

        Matfile = loadmat(path, simplify_cells=True, squeeze_me=True)
        Params = Matfile[Matfile["expname"]]
        uwb_output = uwb_load(Matfile)
        axes = np.squeeze(uwb_output.dta_x)
        data = np.squeeze(uwb_output.dta_ev)
        data = dataset(axes, data, Params)
        data.add_variable(Parameter(name='nAvgs', value=uwb_output.nAvgs))
        return data 


def progress_bar(progress, post=""):

    num_hash = round(progress // 0.05)
    num_space = 20-num_hash

    print(
        "Progress: "+"|"+"#"*num_hash + " " * num_space + "|" + post,
        end='\r')


def progress_bar_frac(num, den):
    
    progress_bar(num/den, f"{num:d} out of {den:d}")
