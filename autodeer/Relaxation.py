from matplotlib.figure import Figure
import matplotlib.cm as cm
import numpy as np
from deerlab import noiselevel
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from autodeer.sequences import Sequence
import deerlab as dl
# ===========================================================================


class CarrPurcellAnalysis():

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

        data = self.data
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
            data = self.data
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
    
    def check_decay(self,level=0.05):
        """
        Checks that the data has decayed to 5% in the time lenght of the dataset.

        Parameters
        ----------
        level : float, optional
            The level to check the decay, by default 0.05
        
        """

        if hasattr(self,"fit_result"):
            decay = self.func(self.axis, *self.fit_result[0]).data
            if decay.min() < level:
                return True
            else:
                return False
        else:
            raise ValueError("No fit result found")

    def find_optimal(
            self, SNR_target, target_time: float, target_step, averages=None) -> float:
        """Calculate the optimal inter pulse delay for a given total measurment
        time. 

        Parameters
        ----------
        SNR_target: float,
            The Signal to Noise ratio target.
        target_time : float
            The target time in hours
        target_shrt : float
            The shot repettition time of target in seconds
        target_step: float
            The target step size in ns.
        averages : int, optional
            The total number of shots taken, by default None. If None, the
            number of shots will be calculated from the dataset.
        

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
        target_step_us = target_step * 1e-3
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

    def calc_optimal_reptime(self, recovery=0.9):
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

cmap = ['#D95B6F','#42A399']

def plot_1Drelax(*args,fig=None, axs=None,cmap=cmap):
    """
    Create a superimposed plot of relaxation data and fits.

    Parameters
    ----------
    args : ad.Analysis
        The 1D relaxation data to be plotted.

    fig : Figure, optional
        The figure to plot to, by default None
    axs : Axes, optional
        The axes to plot to, by default None
    cmap : list, optional
        The color map to use, by default ad.cmap
    
    """

    if fig is None:
        fig, axs = plt.subplots(1,1, figsize=(5,5))
    else:
        axs = fig.subplots(1,1)

    for i,arg in enumerate(args): 
        if arg.dataset.seq_name == 'T2RelaxationSequence':
            xscale = 2
            label='Hahn Echo'
        elif arg.dataset.seq_name == 'CarrPurcellSequence':
            xscale = 4
            label='CP-2'
        elif arg.dataset.seq_name == 'DEERSequence':
            xscale = 4
            label='CP-2'
        else:
            xscale = 4
            label='CP-2'
        
        axs.plot(arg.axis*xscale, arg.data/arg.data.max(), '.', label=label,alpha=0.5,color=cmap[i],mec='none')
        axs.plot(arg.axis*xscale, arg.func(arg.axis,*arg.fit_result[0]), '-',alpha=1,color=cmap[i], lw=2)

    axs.legend()
    axs.set_xlabel('Total Sequence Length / $\mu s$')
    axs.set_ylabel('Signal / $ A.U. $')

    return fig


class RefocusedEcho2DAnalysis():

    def __init__(self, dataset, sequence: Sequence = None) -> None:
        """Analysis and calculation of Refocused Echo 2D data. 

        Parameters
        ----------
        dataset : 
            The dataset to be analyzed.
        sequence : Sequence, optional
            The sequence object describing the experiment. (not currently used)
        """
        self.axis = []
        if 'tau1' in dataset.coords and 'tau2' in dataset.coords:
            self.axis.append(dataset['tau1'])
            self.axis.append(dataset['tau2'])
        elif 'Xt' in dataset.coords:
            self.axis.append(dataset['Xt'])
            self.axis.append(dataset['Yt'])
        
        dataset.epr.correctphasefull
        self.data = dataset.data
        self.dataset = dataset

    def plot2D(self, contour=True, norm = 'Normal', axs=None, fig=None):
        """
        Create a 2D plot of the 2D relaxation data.

        Parameters
        ----------
        contour : bool, optional
            Plot the contour of the data, by default True
        norm : str, optional
            Normalise the data, by default 'Normal'. Options are 'Normal' and 'tau2'. With 'tau2' normalisation, the data is normalised to the maximum of each row.
        axs : Axes, optional
            The axes to plot to, by default None
        fig : Figure, optional 
            The figure to plot to, by default None
        
        """

        data = self.data.real
        if norm == 'Normal':
            data = data/np.max(data)
        elif norm == 'tau2':
            data = data/np.max(data,axis=1)[:,None]

        if axs is None and fig is None:
            fig, axs = plt.subplots()
        elif axs is None:
            axs = fig.subplots(1,1)

        cmap = cm.get_cmap('Purples',lut=None)
        cmap_contour = cm.get_cmap('Greys_r',lut=None)

        axs.pcolormesh(self.axis[0],self.axis[1],data,cmap=cmap)
        if contour is True:
            axs.contour(self.axis[0],self.axis[1],data, cmap=cmap_contour)
        axs.set_xlabel(r'$\tau_1$ / $(\mu s)$')
        axs.set_ylabel(r'$\tau_2$ / $(\mu s)$')
        axs.set_xlim(min(self.axis[0]),max(self.axis[0]))
        axs.set_ylim(min(self.axis[1]),max(self.axis[1]))
        axs.set_aspect('equal')

        return fig

    def plot1D(self,axs=None,fig=None):
        """
        Create a 1D plot of the 2D relaxation data.

        Parameters
        ----------
        axs : Axes, optional
            The axes to plot to, by default None
        fig : Figure, optional
            The figure to plot to, by default None
        """

        # TODO: Expand to include optimal data when the 2D data is not symetrical
        if axs is None and fig is None:
            fig, axs = plt.subplots()
        elif axs is None:
            axs = fig.subplots(1,1)

        data = self.data.real
        data /= np.max(data)

        optimal_4p = np.argmax(data,axis=1)
        
        
        axs.plot(self.axis[0],np.diag(data[:,optimal_4p]),label='4 Pulse',color=cmap[0])
        axs.plot(self.axis[0]*2,np.diag(data),label='5 pulse',color=cmap[1])
        axs.legend()
        axs.set_xlabel(r'$\tau_{evo}$ / $(\mu s)$')
        axs.set_ylabel('Signal / A.U.')

        return fig

        