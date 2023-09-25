import os
import deerlab as dl
import numpy as np
import logging
from autodeer.dataset import Dataset
from autodeer.hardware.ETH_awg_load import uwb_load
from scipy.io import loadmat

log = logging.getLogger('core.Tools')


def eprload(
        path: str, experiment: str = None, type: str = None,
        **kwargs) -> Dataset:
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
    Dataset
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
            return Dataset(t, V)

        else:
            t, V, Params = dl.deerload(
                path, plot=False, full_output=full_output)
            return Dataset(t, V, Params)

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
        # Params = Matfile[Matfile["expname"]]
        if "options" in kwargs:
            uwb_output = uwb_load(Matfile, kwargs["options"])
        else:
            uwb_output = uwb_load(Matfile)
        # axes = uwb_output.dta_x
        # data = uwb_output.dta_ev

        return uwb_output        

def progress_bar(progress, post=""):

    num_hash = round(progress // 0.05)
    num_space = 20-num_hash

    print(
        "Progress: "+"|"+"#"*num_hash + " " * num_space + "|" + post,
        end='\r')


def progress_bar_frac(num, den):
    
    progress_bar(num/den, f"{num:d} out of {den:d}")
