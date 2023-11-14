from autodeer.dataset import Dataset
from autodeer.DEER_analysis import DEERanalysis
import time
import numpy as np
from deerlab.utils import der_snr


class Criteria:
    """
    A class for defining criteria for terminating experiments. This should
    only be subclassed and not used directly.
    """
    def __init__(
            self, name: str, test, description: str = None) -> None:
 
        self.name = name
        self.description = description
        self.test = test
        pass


class TimeCriteria(Criteria):
    def __init__(
            self, name: str, end_time: float, description: str = None,
            ) -> None:
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

        def test_func(Data: Dataset, verbosity=0):
            now = time.time()

            return now > end_time

        super().__init__(name, test_func, description)


class SNRCriteria(Criteria):

    def __init__(
            self, SNR_target: int, description: str = None,verbosity=0) -> None:
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

        def test_func(data: Dataset, verbosity=verbosity):
            # Normalise data
            norm_data = data.data / data.data.max()
            std = der_snr(np.abs(norm_data))
            snr = 1/std
            test = snr > SNR_target
            if verbosity>1:
                print(f"Test: {test}\t - SNR:{snr}")
            return test

        super().__init__("SNR Criteria", test_func, description)


class DEERCriteria(Criteria):

    def __init__(self, mode="Speed", model=None, verbosity=0, update_func=None) -> None:
        """Criteria for running DEER experiments.

        Mode
        ------
        +------------+--------+------+------+-------+
        | Parameter  | Speed  | Low  | Med  | High  |
        +============+========+======+======+=======+
        | MNR        | 20     | 10   | 50   | 100   |
        +------------+--------+------+------+-------+


        Parameters
        ----------
        tau1 : _type_
            _description_
        tau2 : _type_
            _description_
        tau3 : _type_, optional
            _description_, by default None
        mode : str, optional
            _description_, by default "Speed"

        Returns
        -------
        _type_
            _description_
        """
        
        name = "DEERCriteria"
        description = "Criteria for terminating DEER experiments."
        if mode.lower() == "speed":
            MNR_threshold = 20
            regparamrange = (1,1e3)

        elif mode.lower() == "low":
            MNR_threshold = 10
            regparamrange = None
        elif mode.lower() == "med":
            MNR_threshold = 50
            regparamrange = None
        elif mode.lower() == "high":
            MNR_threshold = 100
            regparamrange = None    
        else:
            MNR_threshold = 50
            regparamrange = None

        def test_func(data: Dataset, verbosity=verbosity):
            # fit, _, _ = DEERanalysis(
            #     data.axes[0]/1000 - tau1, data.data,
            #     tau1, tau2, tau3, num_points=100,
            #     compactness=True, precision="Speed", plot=False)
            fit = DEERanalysis(
                data, compactness=False, model=model, regparamrange=regparamrange,verbosity=verbosity,lin_maxiter=50,max_nfev=100
            )
            test = True
            if fit.MNR < MNR_threshold:
                test = False
            
            if update_func is not None:
                update_func(fit)
            
            if verbosity > 0:
                print(f"Test: {test}\t - MNR:{fit.MNR}")
            
            return test
        
        super().__init__(name, test_func, description)
