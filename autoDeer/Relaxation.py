from attr import has
from matplotlib.figure import Figure
import numpy as np
from deerlab import deerload,noiselevel
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

class Carr_Purcell:

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
    
    def fit(self,type = "mono"):

        if type == "mono":
            self.func = lambda x,a,b,e: a*np.exp(-b*x**e)
            p0 = [self.data.max(),1,1]
            bounds = ([0,0,0],[np.inf,np.inf,np.inf])
        else:
            raise ValueError("Type must be one of: mono")
        
        self.popt = curve_fit(self.func,self.axis,self.data,p0=p0,bounds=bounds)
        pass


    def plot(self,norm=True) -> Figure:

        if norm == True:
            data = np.abs(self.data)
            data /= np.max(data)

        fig, ax = plt.subplots()
        if hasattr(self,"popt"):
            ax.plot(self.axis/1000,self.func(self.axis,*self.popt),label='fit')
            ax.plot(self.axis/1000,data,label='data')
            ax.legend()
        else:
            ax.plot(self.axis/1000,data,label='data')


        ax.set_xlabel('Time / us')
        ax.set_ylabel('Normalised Amplitude')
        return fig

    def find_optimal(self,shrt,averages):
        time_per_point = shrt * averages

        data = np.abs(self.data)
        data /= np.max(data)

        self.noise = noiselevel(data)
        data_snr = data/self.noise
        data_snr_avgs = data_snr/np.sqrt(averages)

        # Assume 16ns time step
        dt = 16
        # Target time
        target_time = 2 * 3600

        num_avgs = target_time / (shrt * np.floor(2 * self.axis / 16))

        data_snr_ = data_snr_avgs * np.sqrt(num_avgs)

        self.axis[np.min(np.abs(data_snr_-20))]
