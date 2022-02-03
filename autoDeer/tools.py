import os
import deerlab as dl
from matplotlib.pyplot import plot
import numpy as np
import logging
import h5py

log = logging.getLogger('core.Tools')


def eprload(path:str,experiment:str=None,type=None,**kwargs):

    if type == None: #Use the file ending to guess file type
        filename, file_extension = os.path.splitext(path)

        if file_extension == '.DSC' | '.DTA':
            log.debug('File detected as Bruker')
            type = 'BRUKER'

        elif file_extension == '.h5' | '.hdf5':
            log.debug('File detected as HDF5')
            type = 'HDF5'

        elif file_extension == '.csv' | '.txt':
            log.debug('File detected as csv or txt file')
            type='TXT'

        elif file_extension == '.mat':
            log.debug('File detecetd as Matlab')
            type='MAT'
        
        else:
            log.error("Can't detect file type")
            raise ValueError("Can't detect file type. Please choose the correct file or set type manually \n"+
                            "Valid file types: '.DSC','.DTA','.h5','.hdf5','.csv','.txt','.mat'")
    

    if type == 'BRUKER':
        if 'full_output' in kwargs:
            full_output = kwargs['full_output']
        else:
            full_output = False
        if full_output == False:    
            t,V = dl.deerload(path,plot=False,full_output=full_output)
            return t,V
        else:
            t,V,Params = dl.deerload(path,plot=False,full_output=full_output)
            return t,V,Params

    elif type == 'HDF5':
        f = h5py.File(path,'r')
        groups = list(f.keys())
        
        if experiment == None:
            if 'Spectrometer' in groups:
                groups.remove('Spectrometer')
            exp = f[groups[0]]
        else:
            exp = f[experiment]
        
        attrs = (exp.attrs.keys())
        elements = list(exp.keys())

        if all(item in elements for item in ['Axes','Data']):
            t = exp['Axes']
            V = exp['Data']
        else:
            raise RuntimeError("HDF5 File is missing the datasets ('Axes','Data') in group {experiment}")

        if 'full_output' in kwargs:
            full_output = kwargs['full_output']
        
            if full_output:
                Params = dict(exp.attrs.items())
                return t,V,Params     
            else:
                return t,V
        else:
            return t,V

    elif type == 'TXT':
        if 'full_output' in kwargs:
            full_output = kwargs['full_output']
            del kwargs['full_output']
            if full_output:
                print("WARNING: Can't get metadata from text file")
        data = np.loadtxt(path, *kwargs)
        return data


    elif type == 'MAT':
        pass

