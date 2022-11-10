import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import hilbert
from scipy.io import savemat
import os


# =============================================================================
class dataset:
    """
    The container for all experimental data.
    """
    def __init__(self, axes: np.ndarray, data: np.ndarray, params: dict = None
                 ) -> None:
        """
        Parameters
        ----------
        axes : np.ndarray
            An array of vectors containing the axes of the dataset
        data : np.ndarray
            The data either as an array of n-dimensions, containing either 
            float or complex values.
        params : dict, optional
            A dictionary of experimental parameters, by default None
        """
        self.axes = axes
        self.data = data
        self.params = params
        self.dims = len(self.axes)

        if not np.iscomplexobj(self.data):
            self.data = hilbert(self.data)
        pass

    def plot(self, label: list[str] = None, **kwargs) -> plt.figure:
        """
        Produces a standard quick graph of the data. This is a line plot for 
        1D data and a heatmap for 2D data.

        Parameters
        ----------
        label : list[str], optional
            A list contains labels. [x_label,z_label], by default None

        Returns
        -------
        plt.figure
            _description_
        """
        if self.dims == 1:
            fig, ax = plt.subplots()
            ax.plot(self.axes, np.abs(self.data), label='abs')
            ax.plot(self.axes, np.real(self.data), label='real')
            ax.plot(self.axes, np.imag(self.data), label='imag')
            ax.legend()
            ax.set_ylabel('Signal')
            if label is not None:
                ax.set_xlabel(label[0])
        
    def save(self, file: str) -> None:
        """
        Saves the dataset in a variety of commonly used data formats based of
        the file extension. 

        Extension options: [".mat",".np",".txt"]

        Parameters
        ----------
        file : str
            The file name or path to save to. 
        """
        filename, file_ext = os.path.splittext(file)
        
        if file_ext == ".mat":  # Save as old-style .mat file
            save_dict = {
                "dta": self.data,
                "axes": self.axes,
                "params": self.params
                }
            savemat(filename, save_dict) 
        elif file_ext == ".np":  # Save as numpy file
            if self.dims == 1:
                save_array = np.vstack([self.axes, self.data])
                np.save(filename, save_array)

        elif file_ext == ".txt":  # Save as text file
            if self.dims == 1:
                save_array = np.vstack([self.axes, self.data])
                np.savetxt(filename, save_array, delimiter=",")


# =============================================================================
