from autodeer.classes import Pulse, Parameter
import numpy as np
from scipy.integrate import cumulative_trapezoid


class RectPulse(Pulse):
    """
    Represents a rectangular monochromatic pulse.
    """

    def __init__(
            self, tp, freq, t=None, flipangle=None, **kwargs) -> None:
        Pulse.__init__(
            self, tp=tp, t=t, flipangle=flipangle, **kwargs)
        self.freq = Parameter("freq", freq, "GHz", "Frequency of the Pulse")
        self.Progression = False
        self._buildFMAM(self.func)
        pass

    def func(self, ax):
        nx = ax.shape[0]
        AM = np.ones(nx)
        FM = np.zeros(nx) + self.freq.value
        return AM, FM


# =============================================================================

class HSPulse(Pulse):
    """
    Represents a hyperboilc secant frequency-swept pulse.
    """
    def __init__(self, *, tp, order1, order2, beta, **kwargs) -> None:
        Pulse.__init__(self, tp=tp, **kwargs)
        self.order1 = Parameter(
            "order1", order1, None, "Order 1 of the HS Pulse")
        self.order2 = Parameter(
            "order1", order2, None, "Order 2 of the HS Pulse")
        self.beta = Parameter("beta", beta, None, "Beta of the HS Pulse")

        # Frequency Infomation
        if "BW" in kwargs:
            # Bandwidth + one other
            self.BW = Parameter(
                "Bandwidth", kwargs["BW"], "GHz", "Bandwidth of pulse")
            if "init_freq" in kwargs:
                self.init_freq = Parameter(
                    "init_freq", kwargs["init_freq"], "GHz",
                    "Initial frequency of pulse")
                self.final_freq = Parameter(
                    "final_final", self.BW.value + kwargs["init_freq"], "GHz",
                    "Final frequency of pulse")
            elif "final_freq" in kwargs:
                self.final_freq = Parameter(
                    "final_final", kwargs["final_freq"], "GHz",
                    "Final frequency of pulse")
                self.init_freq = Parameter(
                    "init_freq", kwargs["final_freq"] - self.BW.value, "GHz",
                    "Initial frequency of pulse")
            else:
                raise ValueError()
        
        elif ("init_freq" in kwargs) & ("final_freq" in kwargs):

            self.init_freq = Parameter(
                "f_init", kwargs["init_freq"], "GHz",
                "Initial frequency of pulse")
            self.final_freq = Parameter(
                "f_final", kwargs["final_freq"], "GHz",
                "Final frequency of pulse")
            self.BW = Parameter(
                "Bandwidth", np.abs(kwargs["final_freq"] - kwargs["init_freq"]),
                 "GHz", "Bandwidth of pulse")
        
        else:   
            raise ValueError()
        
        self._buildFMAM(self.func)
        pass


    def func(self, ax):
        beta = self.beta.value
        order1 = self.order1.value
        order2 = self.order2.value
        BW = self.BW.value
        tp = ax.max() - ax.min()
        tcent = tp / 2
        
        nx = ax.shape[0]
        beta_exp1 = np.log(beta*0.5**(1-order1)) / np.log(beta)
        beta_exp2 = np.log(beta*0.5**(1-order2)) / np.log(beta)
        cut = round(nx/2)
        AM = np.zeros(nx)
        AM[0:cut] = 1/np.cosh(
            beta**beta_exp1 * (ax[0:cut]/tp)**order1)
        AM[cut:-1] = 1/np.cosh(
            beta**beta_exp2 * (ax[cut:-1]/tp)**order2)

        FM = BW * cumulative_trapezoid(AM**2,ax,initial=0) /\
             np.trapz(AM**2,ax) + self.init_freq.value

        return AM, FM




# =============================================================================

class ChirpPulse(Pulse):
    """
    Represents a linear frequency-swept pulse.
    """

    def __init__(self, *, tp, **kwargs) -> None:
        Pulse.__init__(self, tp=tp, **kwargs)

        # Frequency Infomation
        if "BW" in kwargs:
            # Bandwidth + one other
            self.BW = Parameter(
                "Bandwidth", kwargs["BW"], "GHz", "Bandwidth of pulse")
            if "init_freq" in kwargs:
                self.init_freq = Parameter(
                    "init_freq", kwargs["init_freq"], "GHz",
                    "Initial frequency of pulse")
            elif "final_freq" in kwargs:
                self.final_freq = Parameter(
                    "final_final", kwargs["final_freq"], "GHz",
                    "Final frequency of pulse")
            else:
                raise ValueError()
        
        elif ("init_freq" in kwargs) & ("final_freq" in kwargs):

            self.init_freq = Parameter(
                "f_init", kwargs["init_freq"], "GHz",
                "Initial frequency of pulse")
            self.final_freq = Parameter(
                "f_final", kwargs["final_freq"], "GHz",
                "Final frequency of pulse")
        
        else:   
            raise ValueError()

        self._buildFMAM(self.func)
        pass

    def func(self, ax):
        nx = ax.shape[0]
        AM = np.ones(nx)

        if hasattr(self, "BW") & hasattr(self, "init_freq"):
            FM = np.linspace(
                self.init_freq.value, self.init_freq.value + self.BW.value, nx)
        elif hasattr(self, "BW") & hasattr(self, "final_freq"):
            FM = np.linspace(
                self.final_freq.value - self.BW.value, self.final_freq.value,
                nx)
        elif hasattr(self, "init_freq") & hasattr(self, "final_freq"):
            FM = np.linspace(
                self.init_freq.value, self.final_freq.value, nx)

        return AM, FM
        

# =============================================================================

class SincPulse(Pulse):
    """
    Represents a sinc shaped monochromatic pulse.
    """

    def __init__(self, *, tp, freq, order, window=None, **kwargs) -> None:
        Pulse.__init__(self, tp=tp, **kwargs)
        self.freq = Parameter("Freq", freq, "GHz", "Frequency of the Pulse")

        self.order = Parameter("Order", order, None, "The sinc pulse order")
        self.window = Parameter(
            "Window", window, None, "The type of window function")

        self._buildFMAM(self.func)
        pass


    def func(self, ax):
        nx = ax.shape[0]
        FM = np.zeros(nx) + self.freq.value
        AM = np.sinc(self.order.value * ax)

        return AM, FM
        