from matplotlib.figure import Figure
import numpy as np
from deerlab import noiselevel
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from autodeer.openepr import dataset

# ===========================================================================


class Carr_Purcell:

    def __init__(self, dataset: dataset) -> None:
        """Analysis and calculation of Carr Purcell decay. 

        Parameters
        ----------
        dataset : dataset
            _description_
        """
        self.axis = dataset.axes
        self.data = dataset.data
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
            self, target_time: float, shrt: float, averages: int) -> float:
        """Calculate the optimal inter pulse delay for a given total measurment
        time. 

        Parameters
        ----------

        target_time : float
            The target s
        shrt : float
            The shot repertition time
        averages : int
            The number of averages in the data, both shots and scans

        Returns
        -------
        float
            The calculated optimal time
        """
        # time_per_point = shrt * averages

        data = np.abs(self.data)
        data /= np.max(data)

        self.noise = noiselevel(data)
        data_snr = data / self.noise
        data_snr_avgs = data_snr / np.sqrt(averages)

        # Assume 16ns time step
        # dt = 16
        # Target time
        target_time = target_time * 3600

        num_avgs = target_time / (shrt * np.floor(2 * self.axis * 1000 / 16))

        data_snr_ = data_snr_avgs * np.sqrt(num_avgs)

        self.optimal = self.axis[np.argmin(np.abs(data_snr_-20))]
        return self.optimal
