from matplotlib.figure import Figure
import numpy as np
from deerlab import noiselevel
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from autodeer.classes import Dataset
from autodeer.sequences import Sequence

# ===========================================================================


class CarrPurcellAnalysis:

    def __init__(self, dataset: Dataset, sequence: Sequence = None) -> None:
        """Analysis and calculation of Carr Purcell decay. 

        Parameters
        ----------
        dataset : Dataset
            _description_
        """
        self.axis = dataset.axes[0]
        self.data = dataset.data
        self.seq = sequence
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
            # bounds = ([0, 0, -10, 0],[2, 1000, 10, 1000])
        else:
            raise ValueError("Type must be one of: mono")
        
        self.fit_result = curve_fit(self.func, self.axis, data, p0=p0)
        return self.fit_result

    def plot(self, norm: bool = True) -> Figure:
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

        fig, ax = plt.subplots()
        if hasattr(self, "fit_result"):
            ax.plot(self.axis, data, '.', label='data', color='0.6', ms=6)
            ax.plot(self.axis, self.func(
                self.axis, *self.fit_result[0]), label='fit', color='C1', lw=2)

            ax.legend()
        else:
            ax.plot(self.axis, data, label='data')

        ax.set_xlabel('Time / us')
        ax.set_ylabel('Normalised Amplitude')
        return fig

    def find_optimal(
            self, averages: int, SNR_target, target_time: float, target_shrt: float, target_step) -> float:
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
        data = np.abs(self.data)
        data /= np.max(data)
        if hasattr(self,"fit_result"):
            calc_data = self.func(self.axis,*self.fit_result[0])
        else:
            calc_data = data

        # averages = self.seq.shots.value * self.seq.averages.value
        self.noise = noiselevel(data)
        data_snr = calc_data / self.noise
        data_snr_avgs = data_snr / np.sqrt(averages)

        # Target time
        target_time = target_time * 3600

        g = (target_time * target_step / target_shrt) * 1/(self.axis)
        f = (SNR_target/data_snr_avgs)**2

        self.optimal = self.axis[np.argmin(np.abs(g-f))]
        return self.optimal
