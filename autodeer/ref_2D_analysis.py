from pyepr import Sequence
import numpy as np
from deerlab import noiselevel
from scipy.linalg import svd
import matplotlib.pyplot as plt
from autodeer.colors import primary_colors

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
        
        dataset = dataset.epr.correctphasefull
        self.data = dataset.data
        self.dataset = dataset
        
        data = self.data / np.max(self.data)
        self.noise = noiselevel(data)
        

    def _smooth(self,elements=3):
        """
        Used SVD to smooth the 2D data.

        Parameters
        ----------
        elements : int, optional
            The number of elements to use in the smoothing, by default 3
        
        Returns
        -------
        np.ndarray
            The smoothed data.
        """

        U,E,V = svd(self.data.real)
        E_mod = E.copy()
        E_mod[elements:] = 0
        mod_data = U @ np.diag(E_mod) @ V
        self.data_smooth = mod_data

        return self.data_smooth
    
    def _calc_optimal_4p_values(self):
        if not hasattr(self,'data_smooth'):
                self._smooth()
        data = self.data_smooth
        data /= np.max(data)

        optimal_4p = np.argmax(data,axis=0)
        self.optimal_4p_V = np.max(data,axis=0)
        self.optimal_4p_tau1 = self.axis[0].data[optimal_4p]
        self.optimal_4p_tau2 = self.axis[1].data
        

    def plot2D(self, contour=True,smooth=False, norm = 'Normal', axs=None, fig=None):
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
        if smooth is True:
            if not hasattr(self,'data_smooth'):
                self._smooth()
            data = self.data_smooth
        else:
            data = self.data.real

        if norm == 'Normal':
            data = data/np.max(data)
        elif norm == 'tau2':
            data = data/np.max(data,axis=1)[:,None]

        if axs is None and fig is None:
            fig, axs = plt.subplots()
        elif axs is None:
            axs = fig.subplots(1,1)

        cmap = plt.get_cmap('Purples',lut=None)
        cmap_contour = plt.get_cmap('Greys_r',lut=None)


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

        if not hasattr(self,'data_smooth'):
                self._smooth()
        data = self.data_smooth
        data /= np.max(data)

        optimal_4p = np.argmax(data,axis=1)
        
        
        axs.plot(self.axis[0],np.diag(data[:,optimal_4p]),label='4 Pulse',color=primary_colors[0])
        axs.plot(self.axis[0]*2,np.diag(data),label='5 pulse',color=primary_colors[1])
        axs.legend()
        axs.set_xlabel(r'$\tau_{evo}$ / $(\mu s)$')
        axs.set_ylabel('Signal / A.U.')

        return fig
    
    def find_optimal(self,type:str,SNR_target, target_time: float, target_step, averages=None) -> float:
        """Calculate the optimal inter pulse delay for a given total measurment time, using either 4pulse or 5pulse data.

        Parameters
        ----------
        type : str
            The type of data to use, either '4pDEER' or '5pDEER'
        SNR_target : float
            The Signal to Noise ratio target.
        target_time : float
            The target time in hours
        target_step: float
            The target step size in ns.
        averages : int, optional
            The total number of shots taken, by default None. If None, the
            number of shots will be calculated from the dataset.

        Returns
        -------
        tau1: float
            The calculated optimal tau1 in us
        tau2: float
            The calculated optimal tau2 in us

        Notes:
        ------
        The shot repetition time is assumed to be the same as for the relxaation data and is taken from the dataset.
        """

        dataset = self.dataset
        if averages is None:
            averages = dataset.nAvgs * dataset.shots * dataset.nPcyc
        target_shrt = dataset.reptime * 1e-6

        if not hasattr(self,'data_smooth'):
                self._smooth()
        data = self.data_smooth
        raw_data = np.abs(self.data)
        raw_data /= np.max(raw_data)
        data /= np.max(data)

        calc_data = data

        # averages = self.seq.shots.value * self.seq.averages.value
        self.noise = noiselevel(raw_data)
        data_snr = calc_data / self.noise
        data_snr_avgs = data_snr / np.sqrt(averages)

        if type == '4pDEER':
            data_snr_avgs_tau2 = np.max(data_snr_avgs,axis=1)

            # Target time
            target_time = target_time * 3600
            g = (target_time * target_step / target_shrt) * 1/(self.axis[1].data)
            f = (SNR_target/data_snr_avgs_tau2)**2

            tau2_idx = np.argmin(np.abs(g-f))
            self.tau2 = self.axis[1].data[tau2_idx]
            self.tau1 = self.axis[0].data[np.argmax(data_snr_avgs[:,tau2_idx])]
            return self.tau1, self.tau2
        
        elif type == '5pDEER':
            data_snr_avgs_CP = np.diag(data_snr_avgs)
            target_time = target_time * 3600
            g = (target_time * target_step / target_shrt) * 1/(self.axis[1].data)
            f = (SNR_target/data_snr_avgs_CP)**2
            tau2_idx = np.argmin(np.abs(g-f))
            self.tau2 = self.axis[1].data[tau2_idx]
            return self.tau2, self.tau2
        

    
    def optimal_tau1(self,tau2=None,):
        if not hasattr(self,'data_smooth'):
                self._smooth()
        data = self.data_smooth
        tau2_idx = np.argmin(np.abs(self.axis[1].data - tau2))
        self.tau1 = self.axis[0].data[np.argmax(data[:,tau2_idx])]
        return self.tau1          

    def __call__(self, x, y, norm=True, SNR=False):
        """
        Evaluate the fit or data at a given x value.

        Parameters
        ----------
        x : float
            The x value to evaluate the data at.
        norm : bool, optional
            Normalise the data to the maximum, by default True
        SNR : bool, optional
            Return the SNR_per_sqrt(shot) for this data point, by default False
        source : str, optional
            The source of the data, either 'fit' or 'data', by default None
            If None, the source is determined by the presence of a fit result.
        
        """
        
        x_idx = np.abs(self.axis[0] - x).argmin()
        y_idx = np.abs(self.axis[1] - y).argmin()
        V = self.data[y_idx,x_idx]

        if norm is True:
            V /= np.max(self.data)

        if SNR is True:
            V /= self.noise
            V /= np.sqrt(self.dataset.nAvgs * self.dataset.shots * self.dataset.nPcyc)

        return V

