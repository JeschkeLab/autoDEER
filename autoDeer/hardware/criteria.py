from autoDeer.hardware.openepr import dataset
import time
import numpy as np
from deerlab.utils import der_snr


class Criteria:

    def __init__(
            self, name: str, test, description: str = None) -> None:
        
        self.name = name
        self.description = description
        self.test = test
        pass


class TimeCriteria(Criteria):
    def __init__(
            self, name: str, end_time: float, description: str = None) -> None:
        """Criteria testing for a specific finishing time. The finishing time 
        is given as absolute time in the locale of the computer, it is *not* 
        how the long the measurment continues for. 

        Parameters
        ----------
        name : str
            _description_
        end_time : float
            Finishing time
        description : str, optional
            _description_, by default None
        """

        def test_func(Data: dataset):
            now = time.time()

            return now > end_time
                
        super().__init__(name, test_func, description)


class SNRCriteria(Criteria):

    def __init__(
            self, name: str, SNR_target: int, description: str = None) -> None:
        """Criteria testing for signal to noise ratio. This checks the SNR of 
        the normalised absolute data using the deerlab SNR noise estimation
        which is based on the work by Stoher et. al. [1]

        Parameters
        ----------
        name : str
            _description_
        SNR_target : int
            The mimimum SNR value.
        description : str, optional
            _description_, by default None

        References
        -----------
        [1] F. Stoehr, R. White, M. Smith, I. Kamp, R. Thompson, D. Durand,
        W. Freudling, D. Fraquelli, J. Haase, R. Hook, T. Kimball, M. Kummel,
        K. Levay, M. Lombardi, A. Micol, T. Rogers DERSNR: A Simple & General
        Spectroscopic Signal-to-Noise Measurement Algorithm Astronomical Data
        Analysis Software and Systems XVII, ASP Conference Series, Vol. 30,
        2008, p5.4
        """

        def test_func(data: dataset):
            # Normalise data
            norm_data = data.data / data.data.max()
            std = der_snr(np.abs(norm_data))
            snr = 1/std
            return snr > SNR_target

        super().__init__(name, test_func, description)
