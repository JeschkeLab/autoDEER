from autodeer.DEER_analysis import DEERanalysis
import time
import numpy as np
from deerlab.utils import der_snr
from deerlab import noiselevel
import logging

log = logging.getLogger('autoDEER.criteria')


class Criteria:
    """
    A class for defining criteria for terminating experiments. This should
    only be subclassed and not used directly.
    """
    def __init__(
            self, name: str, test, description: str = '', end_signal=None) -> None:
 
        self.name = name
        self.description = description
        self.test = test
        self.end_signal = end_signal
        pass

    def __add__(self, __o:object):

        if not isinstance(__o,Criteria):
            raise RuntimeError("Only objects of the `Criteria` class can be summed together")
        
        new_name = self.name +' + ' + __o.name
        new_desc = self.description + ' + ' + __o.description

        def new_func(data, verbosity=0):
            test1 = self.test(data,verbosity)
            test2 = __o.test(data,verbosity)
            test = test1 or test2
            test_msg = f"Test {new_name}: {test}"
            log.debug(test_msg)
            return test
        
        if callable(self.end_signal) and callable(__o.end_signal):
            def end_signal():
                self.end_signal()
                __o.end_signal()
        elif callable(self.end_signal):
            end_signal = self.end_signal
        elif callable(__o.end_signal):
            end_signal = __o.end_signal
        else:
            end_signal = None
        
        new_crit = Criteria(new_name,new_func,new_desc,end_signal=end_signal)

        return new_crit


class TimeCriteria(Criteria):
    def __init__(
            self, name: str, end_time: float, description: str = '',
            **kwargs) -> None:
        """Criteria testing for a specific finishing time. The finishing time 
        is given as absolute time in the locale of the computer, it is *not* 
        how the long the measurment continues for. 

        Parameters
        ----------
        name : str
            _description_
        end_time : float
            Finishing time in seconds since epoch
        description : str, optional
            _description_, by default None
        """

        def test_func(Data, verbosity=0):
            now = time.time()

            return now > end_time

        super().__init__(name, test_func, description,**kwargs)


class SNRCriteria(Criteria):

    def __init__(
            self, SNR_target: int, description: str = '',verbosity=0,**kwargs) -> None:
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

        def test_func(data, verbosity=verbosity):
            # Normalise data
            norm_data = data.data / data.data.max()
            std = der_snr(np.abs(norm_data))
            snr = 1/std
            test = snr > SNR_target
            test_msg = f"Test {self.name}: {test}\t - SNR:{snr}"
            log.debug(test_msg)
            if verbosity>1:
                print(test_msg)
            return test

        super().__init__("SNR Criteria", test_func, description,**kwargs)


class DEERCriteria(Criteria):

    def __init__(self, mode="Speed", model=None, verbosity=0, update_func=None,**kwargs) -> None:
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

        def test_func(data, verbosity=verbosity):
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
            test_msg = f"Test {self.name}: {test}\t - MNR:{fit.MNR}"
            log.debug(test_msg)
            if verbosity > 0:
                print(test_msg)
            
            return test
        
        super().__init__(name, test_func, description,**kwargs)
