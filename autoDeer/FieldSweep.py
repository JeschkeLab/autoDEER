from matplotlib.figure import Figure
import numpy as np
import matplotlib.pyplot as plt
from autodeer.hardware.openepr import dataset


class FieldSweep():

    def __init__(self, dataset: dataset) -> None:
        """Analysis and calculation of FieldSweep Experiment. 

        Parameters
        ----------
        dataset : dataset
            _description_
        """
        self.axis = dataset.axes
        self.data = dataset.data
        pass

    def find_max(self) -> float:
        """Calculates the maximum field

        Returns
        -------
        float
            Max field
        """
        max_id = np.argmax(np.abs(self.data))
        self.max_field = self.axis[max_id]

        return self.max_field

    def calc_gyro(self, det_frq: float) -> float:
        """Calculates the gyromagnetic ratio for a given frequency

        Parameters
        ----------
        det_frq : float
            The detection frequency for the field sweep.

        Returns
        -------
        float
            The gyromagnetic ratio in G/GHz.
        """

        if not hasattr(self, "max_field"):
            self.find_max()
        self.det_frq = det_frq
        self.gyro = det_frq/self.max_field
        hf_x = det_frq - self.gyro*self.axis
        self.fs_x = det_frq + hf_x
        self.fs_x = det_frq - self.gyro*self.axis
        return self.gyro

    def plot(self, norm: bool = True, axis: str = "time") -> Figure:
        """Generate a field sweep plot

        Parameters
        ----------
        norm : bool, optional
            Nomarlisation of the plot to a maximum of 1, by default True
        axis : str, optional
            plot field sweep on either the "time" axis or "freq" axis

        Returns
        -------
        plt.figure
            matplotlib figure
        """

        if norm is True:
            data = self.data
            data /= np.max(np.abs(data))

        if axis.lower() == 'time':
            fig, ax = plt.subplots()
            ax.plot(self.axis, np.abs(data), label='abs')
            ax.plot(self.axis, np.real(data), label='real')
            ax.plot(self.axis, np.imag(data), label='imag')
            ax.legend()
            ax.set_xlabel('Field G')
            ax.set_ylabel('Normalised Amplitude')
            if hasattr(self, "max_field"):
                min_value = np.min(np.hstack([np.real(data), np.imag(data)]))
                max_value = np.min(np.hstack([np.real(data), np.imag(data)]))
                ax.vlines(
                    self.max_field, min_value, max_value, label="Maximum")
        elif axis.lower() == 'freq':
            if not hasattr(self, "fs_x"):
                raise RuntimeError("Please run fieldsweep.calc_gyro() first")
            fig, ax = plt.subplots()
            ax.plot(self.fs_x, np.abs(data), label='abs')
            ax.set_xlabel('Frequency GHz')
            ax.set_ylabel('Normalised Amplitude')

        return fig
