from autodeer.DEER_analysis import DEERanalysis
import logging
from pyepr.criteria import Criteria

log = logging.getLogger('autoDEER.criteria')

class DEERCriteria(Criteria):

    def __init__(self, mode="Speed", model=None, verbosity=0, update_func=None,**kwargs) -> None:
        """Criteria for running DEER experiments.

        Mode
        ------
        +------------+--------+------+------+-------+
        | Parameter  | Speed  | Low  | Med  | High  |
        +============+========+======+======+=======+
        | MNR        | 20     | 10   | 50   | 150   |
        +------------+--------+------+------+-------+


        Parameters
        ----------
        
        mode : str, optional
            _description_, by default "Speed"

        Returns
        -------
        _type_
            _description_
        """
        
        name = "DEERCriteria"
        description = "Criteria for terminating DEER experiments."
        if isinstance(mode, (int,float)):
            MNR_threshold = mode
            regparamrange = None
        elif mode.lower() == "speed":
            MNR_threshold = 20
            regparamrange = (1,1e3)

        elif mode.lower() == "low":
            MNR_threshold = 10
            regparamrange = None
        elif mode.lower() == "med":
            MNR_threshold = 50
            regparamrange = None
        elif mode.lower() == "high":
            MNR_threshold = 150
            regparamrange = None    
        else:
            MNR_threshold = 50
            regparamrange = None

        if 'compactness' in kwargs:
            compactness = kwargs.pop('compactness')
        else:
            compactness = False

        def test_func(data, verbosity=verbosity):
            # fit, _, _ = DEERanalysis(
            #     data.axes[0]/1000 - tau1, data.data,
            #     tau1, tau2, tau3, num_points=100,
            #     compactness=True, precision="Speed", plot=False)



            fit = DEERanalysis(
                data, compactness=True, model=model, regparamrange=regparamrange,verbosity=verbosity,lin_maxiter=50,max_nfev=100
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


# class BackgroundCriteria(Criteria):

#     def __init__(self, SNR_target:float, model=None, verbosity=0, update_func=None, **kwargs)->None:
#         """Criteria for terminating DEER background experiments.

#         Parameters
#         ----------
#         SNR_target : float
#             Target signal-to-noise ratio.

#         Returns
#         -------
#         _type_
#             _description_
#         """
#         name = "BackgroundCriteria"
#         description = "Criteria for terminating background experiments."

#         def test_func(data, verbosity=verbosity):

#             fit = BackgroundAnalysis(data, model=model, verbosity=verbosity,lin_maxiter=50,max_nfev=100)
            
#             if fit.SNR < SNR_target:
#                 test = False
            
#             if update_func is not None:
#                 update_func(fit)

#             test_msg = f"Test {self.name}: {test}\t - SNR:{fit.SNR}"
#             log.debug(test_msg)
#             if verbosity > 0:
#                 print(test_msg)
            
#             return test



#         super().__init__(name, test_func, description,**kwargs)