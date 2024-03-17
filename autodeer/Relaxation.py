from matplotlib.figure import Figure
import numpy as np
from deerlab import noiselevel
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from autodeer.sequences import Sequence
import deerlab as dl
# ===========================================================================


class CarrPurcellAnalysis:

    def __init__(self, dataset, sequence: Sequence = None) -> None:
        """Analysis and calculation of Carr Purcell decay. 

        Parameters
        ----------
        dataset : 
            _description_
        """
        # self.axis = dataset.axes[0]
        # self.data = dataset.data
        if 'tau1' in dataset.coords:
            self.axis = dataset['tau1']
        elif 'tau' in dataset.coords:
            self.axis = dataset['tau']
        elif 't' in dataset.coords:
            self.axis = dataset['t']
        elif 'step' in dataset.coords:
            self.axis = dataset['step'] 
        else:
            self.axis = dataset['X']
        
        dataset.epr.correctphasefull
        self.data = dataset.data
        self.dataset = dataset
        # if sequence is None and hasattr(dataset,'sequence'):
        #     self.seq = dataset.sequence
        # else:
        #     self.seq = sequence
        pass
    
    def fit(self, type: str = "mono"):
        """Fit the experimental CP decay

        Parameters
        ----------
        type : str, optional
            Either a mono or double exponential decay model, by default "mono"

        """

        data = np.abs(self.data)
        data /= np.max(data)

        if type == "mono":
            self.func = lambda x, a, b, e: a*np.exp(-b*x**e)
            p0 = [1, 1, 2]
            bounds = ([0, 0, 0],[2, 1000, 10])
        elif type == "double":
            self.func = lambda x, a, b, e, c, d, f: a*np.exp(-b*x**e) + c*np.exp(-d*x**f)
            p0 = [1, 1, 2, 1, 1, 2]
            bounds = ([0, 0, 0, 0, 0, 0],[2, 1000, 10, 2, 1000, 10])
        else:
            raise ValueError("Type must be one of: mono")
        
        self.fit_type = type
        self.fit_result = curve_fit(self.func, self.axis, data, p0=p0, bounds=bounds)
        return self.fit_result

    def plot(self, norm: bool = True, axs=None, fig=None) -> Figure:
        """Plot the carr purcell decay with fit, if avaliable.

        Parameters
        ----------
        norm : bool, optional
            Normalise the fit to a maximum of 1, by default True
        

        Returns
        -------
        Figure
            The figure.
        """

        if norm is True:
            data = np.abs(self.data)
            data /= np.max(data)

        if axs is None and fig is None:
            fig, axs = plt.subplots()

        if hasattr(self, "fit_result"):
            axs.plot(self.axis, data, '.', label='data', color='0.6', ms=6)
            axs.plot(self.axis, self.func(
                self.axis, *self.fit_result[0]), label='fit', color='C1', lw=2)

            axs.legend()
        else:
            axs.plot(self.axis, data, label='data')

        axs.set_xlabel('Time / us')
        axs.set_ylabel('Normalised Amplitude')
        return fig

    def find_optimal(
            self, SNR_target, target_time: float, target_step, averages=None) -> float:
        """Calculate the optimal inter pulse delay for a given total measurment
        time. 

        Parameters
        ----------

        averages : int
            The number of averages in the data, shots, scans and phase cycles
        SNR_target: float,
            The Signal to Noise ratio target.
        target_time : float
            The target time in hours
        target_shrt : float
            The shot repettition time of target in seconds
        target_step: float
            The target step size in us.
        

        Returns
        -------
        float
            The calculated optimal time in us
        """
        # time_per_point = shrt * averages
        dataset = self.dataset
        if averages is None:
            averages = dataset.nAvgs * dataset.shots * dataset.nPcyc
        target_shrt = dataset.reptime * 1e-6

        data = np.abs(self.data)
        data /= np.max(data)

        if hasattr(self,"fit_result"):
            calc_data = self.func(self.axis.data,*self.fit_result[0])
        else:
            calc_data = data

        # averages = self.seq.shots.value * self.seq.averages.value
        self.noise = noiselevel(data)
        data_snr = calc_data / self.noise
        data_snr_avgs = data_snr / np.sqrt(averages)

        # Target time
        target_time = target_time * 3600

        g = (target_time * target_step / target_shrt) * 1/(self.axis.data)
        f = (SNR_target/data_snr_avgs)**2

        self.optimal = self.axis.data[np.argmin(np.abs(g-f))]
        return self.optimal

class ReptimeAnalysis():

    def __init__(self, dataset, sequence: Sequence = None) -> None:
        """Analysis and calculation of Reptime based saturation recovery. 

        Parameters
        ----------
        dataset :
            The dataset to be analyzed.
        sequence : Sequence, optional
            The sequence object describing the experiment. (not currently used)
        """
        # self.axis = dataset.axes[0]
        self.axis = dataset['reptime']
        # if self.axis.max() > 1e4:
        #     self.axis /= 1e3 # ns -> us
        # self.data = dataset.data/np.max(dataset.data)
        
        if np.iscomplexobj(dataset.data):
            self.data = dataset.epr.correctphase
        else:
            self.data = dataset

        self.data.data /= np.max(self.data.data)
        self.seq = sequence
        pass

    def fit(self, **kwargs):
        def func(t,A,T1):
            return A*(1-np.exp(-t/T1))
        self.func = func

        # mymodel = dl.Model(func,constants='t')
        # mymodel.T1.set(lb=0,ub=np.inf,par0=1.8e3)
        # mymodel.T1.unit = 'us'
        # mymodel.T1.description = 'T1 time'

        # results = dl.fit(mymodel,self.data.real,self.axis,reg=False,**kwargs)
        # self.fit_result = results

        self.fit_result = curve_fit(func, self.axis, self.data, p0=[1,1.8e3])

        return self.fit_result

    def plot(self, axs=None, fig=None):

        if axs is None and fig is None:
            fig, axs = plt.subplots()

        axs.plot(self.axis, self.data, '.', label='data', color='0.6', ms=6)

        if hasattr(self,'fit_result'):
            axs.plot(self.axis, self.func(self.axis,*self.fit_result[0]), label='fit', color='C1', lw=2)
            axs.vlines(self.fit_result[0][1],0,1,linestyles='dashed',label='T1 = {:.3g} us'.format(self.fit_result[0][1]),colors='C0')

            if hasattr(self,'optimal'):
                axs.vlines(self.optimal,0,1,linestyles='dashed',label='optimal = {:.3g} us'.format(self.optimal),colors='C2')

        axs.set_xlabel('reptime (us)')
        axs.set_ylabel('normalized signal')
        axs.legend()
        return fig

    def calc_optimal_reptime(self, recovery=0.8):
        # Calculates the x% recovery time
        self.optimal = self.fit_result[0][1]*np.log(1/(1-recovery))
        return self.optimal

def detect_ESEEM(dataset,type='deuteron', threshold=1.5):
    """Detect if the dataset is an ESEEM experiment.

    Parameters
    ----------
    dataset : xr.DataArray
        The dataset to be analyzed.
    
    type : str, optional
        The type of ESEEM experiment, either deuteron or proton, by default 'deuteron'
    
    threshold : float, optional
        The SNR threshold for detection, by default 1.5

    Returns
    -------
    bool
        True if ESEEM is detected, False if not.
    """
    

    D_freq = 4.10663 * dataset.B *1e-4 *np.pi /2
    P_freq = 26.75221 * dataset.B *1e-4 *np.pi /2

    def find_pnl(freq):
        fft_data = np.abs(dataset.epr.fft)
        index = np.abs(fft_data.X - freq).argmin().data

        peak = 2 /fft_data.size * fft_data[index]

        noiselevel = 2/fft_data.size * fft_data[index-8:index+8].mean()

        return peak/noiselevel
        
    if type == 'deuteron':
        peak = find_pnl(D_freq)
    elif type == 'proton':
        peak = find_pnl(P_freq)
    else:
        raise ValueError('type must be deuteron or proton')

    if peak > threshold:
        return True
    else:
        return False

    
    
