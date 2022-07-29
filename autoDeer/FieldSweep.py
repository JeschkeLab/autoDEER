from matplotlib.figure import Figure
import numpy as np
from deerlab import deerload
import matplotlib.pyplot as plt


class FieldSweep():

    def __init__(self) -> None:
        pass

    def import_from_bruker(self,file_path)->None:

        t,V = deerload(file_path, False, True)
        self.axis = t
        self.data = V
        pass 

    def import_from_numpy(self,axis:np.ndarray,spectrum:np.ndarray) -> None:
        self.axis = axis
        self.data = spectrum
        pass

    def import_from_dataclass(self,dataclass) -> None:
        self.axis = dataclass.time
        self.data = dataclass.data
        pass

    def find_max(self) -> float:
        max_id = np.argmax(np.abs(self.data))
        self.max_field = self.axis[max_id]

        return self.max_field

    def calc_gyro(self,det_frq):

        if not hasattr(self,"max_field"):
            self.find_max()
        
        self.gyro = det_frq/self.max_field
        hf_x = det_frq - self.gyro*self.axis
        self.fs_x = det_frq + hf_x
        self.fs_x = det_frq - self.gyro*self.axis
        return self.gyro

    def plot(self,norm=True)->Figure:

        if norm == True:
            data = self.data
            data /= np.max(np.abs(data))


        fig, ax = plt.subplots()
        ax.plot(self.axis,np.real(data),label='real')
        ax.plot(self.axis,np.imag(data),label='imag')
        ax.legend()
        ax.set_xlabel('Field G')
        ax.set_ylabel('Normalised Amplitude')
        if hasattr(self,"max_field"):
            min_value = np.min(np.hstack([np.real(data),np.imag(data)]))
            max_value = np.min(np.hstack([np.real(data),np.imag(data)]))
            ax.vlines(self.max_field,min_value,max_value,label="Maximum")

        return fig