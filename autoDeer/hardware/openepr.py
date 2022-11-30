import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import hilbert
from scipy.io import savemat
import os
from autoDeer import __version__


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
# Super Classes Pulses
# =============================================================================

class Sequence:

    def __init__(self, **kwargs) -> None:

        self.pulses = []
        self.num_pulses = len(self.pulses)

        if "B" in kwargs:
            self.B = Parameter(
                "Field", kwargs["B"], "Gauss",
                "The static B0 for the experiment")
        
        if "LO" in kwargs:
            self.LO = Parameter(
                "LO", kwargs["LO"], "GHz",
                "The local oscillator frequency.")
        
        if "reptime" in kwargs:
            self.reptime = Parameter(
                "reptime", kwargs["reptime"], "us",
                "The shot repetition time")
        
        if "averages" in kwargs:
            self.reptime = Parameter(
                "Averages", kwargs["averages"], "None",
                "The number of averages to perform.")
        
        if "shots" in kwargs:
            self.reptime = Parameter(
                "Shots", kwargs["shots"], "None",
                "The number of shots per scan.")
        pass

    def plot(self) -> None:

        pass

    def addPulse(self, pulse):
        if type(pulse) == Pulse or issubclass(type(pulse), Pulse):
            self.pulses.append(pulse)
        
        elif type(pulse) == list:
            for el in pulse:
                self.pulses.append(el)
        self.num_pulses = len(self.pulses)

    def _buildPhaseCycle(self):
        # Identify pulses which are phase cycled

        pcyc_pulses = []
        pulse_cycles = []
        det_cycles = []
        for ix, pulse in enumerate(self.pulses):
            if pulse.pcyc is not None:
                pcyc_pulses.append(ix)
                pulse_cycles.append(np.array(pulse.pcyc[0]))
                det_cycles.append(np.array(pulse.pcyc[1]))
        
        # Build expanded phase cycle
        func = lambda x: np.arange(0, len(x))
        map(func, pulse_cycles)
        n_pulses = len(pulse_cycles)

        m = list(map(func, pulse_cycles))
        grids = np.meshgrid(*m, indexing='ij')
        expanded_cycles = []
        expanded_dets = []

        for i in range(0, n_pulses):
            expanded_cycles.append(
                pulse_cycles[i][grids[i].flatten(order='F')])
            expanded_dets.append(det_cycles[i][grids[i].flatten(order='F')])

        self.pcyc_vars = pcyc_pulses
        self.pcyc_cycles = np.stack(expanded_cycles)
        self.pcyc_dets = np.prod(np.stack(expanded_dets), axis=0)
        self.pcyc_n = self.pcyc_cycles.shape[1]
        return self.pcyc_vars, self.pcyc_cycles, self.pcyc_dets

    def _buildProgTable(self):
        #           id, pulse #, variable, axis
        progTable = [[], [], [], []]
        for n, pulse in enumerate(self.pulses):
            # Loop over all pulses
            for var_name in vars(pulse):
                var = getattr(pulse, var_name)
                # Loop over all pulse parameters
                if type(var) is Parameter:
                    if var.progressive is True:
                        for i in range(0, len(var.prog)):
                            progTable[0].append(var.prog[i][0]) 
                            progTable[1].append(n)
                            progTable[2].append(var_name) 
                            progTable[3].append(var.prog[i][1]) 

        # TODO: Add support for non_pulse parameters
        #             
        self.progTable = progTable
        return self.progTable
        

    def addPulseProg(self, pulse_id, variable, axis_id, axis) -> None:
        var = getattr(self.pulses[pulse_id], variable)
        var.add_progression(axis_id, axis)
        pass
    
    def addPulsesProg(
            self, pulses, variables, axis_id, axis, multipliers=None) -> None:
        
        if multipliers is None:
            multipliers = np.ones(len(pulses))

        for i, pulse in enumerate(pulses):
            self.addPulseProg(pulse, variables[i], axis_id,
                              axis * multipliers[i])                  

    def _checkRect(self) -> bool:
        """Checks if all the pulses in the sequence are rectangular.
        """
        test = True

        for pulse in self.pulses:
            if type(pulse) is not RectPulse:
                test = False

        return test


# =============================================================================


class Pulse:

    def __init__(self, t, tp, scale) -> None:
        self.t = Parameter("Time", t, "ns", "Start time of pulse")
        self.tp = Parameter("Length", tp, "ns", "Length of the pulse")
        self.scale = Parameter("Scale", scale, None, "Amplitude of pulse")
        self.Progression = False
        self.pcyc = None
        pass

    def _addPhaseCycle(self, phases, detections=None):
        if detections is None:
            detections = np.ones(len(phases))
        self.pcyc = [phases, detections]
        pass

    def _buildFMAM(self, func):
        self.ax = np.linspace(-self.tp.value/2, self.tp.value/2, 1000)
        self.AM, self.FM = func(self.ax)
        self.complex = self.AM *\
            (np.cos(self.FM*self.ax) + 1j * np.sin(self.FM*self.ax))
        self.FM
        self.AM
        return self.AM, self.FM

    def plot(self, pad=1000):
        dt = self.ax[1] - self.ax[0]
        tax = self.t.value + self.ax
        fwd_ex = np.linspace(tax[0] - dt * pad, tax[0], pad)
        rev_ex = np.linspace(
            tax[-1] + dt, tax[-1] + dt * pad, pad, endpoint=True)
        tax = np.concatenate((fwd_ex, tax, rev_ex))
        data = np.pad(self.complex, pad)

        plt.plot(tax, np.real(data), label="Re")
        plt.plot(tax, np.imag(data), label="Im")
        plt.legend()
        plt.xlabel("Time (ns)")
        plt.ylabel("Amplitude (A.U.)")
        pass

    def plot_fft(self):
        pass

    def __str__(self):
        # Build header line
        header = "#" * 79 + "\n" + "OpenEPR Pulse Definition" + \
                 "\n" + "#" * 79 + "\n"
        
        # Build Overviews
            
        overview_params = "Time Pos (ns) \tLength (ns) \tScale (0-1)" +\
            "\tProgressive \n"
        
        if type(self) is Detection:
            overview_params += f"{self.t.value}" + "\t\t" + \
                f"{self.tp.value}" + "\t\t" + "\t\t" +\
                f"{self.Progression}" + "\n" + "#" * 79 + "\n"        
        else:
            overview_params += f"{self.t.value}" + "\t\t" + \
                f"{self.tp.value}" + "\t\t" + f"{self.scale.value}" + "\t\t" +\
                f"{self.Progression}" + "\n" + "#" * 79 + "\n"

        # Build Pulse Parameter Table
        param_table = "Name \t\t Value \t Unit \t Description\n"
        param_table += "---- \t\t ----- \t ---- \t -----------\n"
        for attr_name in dir(self):
            attr = getattr(self, attr_name)
            if type(attr) is Parameter:
                param_table += f"{attr.name} \t\t {attr.value} \t" +\
                    f"{attr.unit} \t {attr.description} \n"

        # Build Pulse Progression Table
        if self.Progression is True:
            prog_table = "#" * 79 + "\n" + "Progressive Variables:" + "\n"
            prog_table += "Variable \t Axis \t Length \n"
            prog_table += "-------- \t ---- \t ------ \n"
            prog_table += f"{self.Prog_var} \t\t {self.Prog_axis} \t" +\
                f"{self.Prog_n} \n"
        else:
            prog_table = ""

        # Build Pulse Phase Cycle Table
        if self.pcyc is not None:
            pcyc_table = "#" * 79 + "\n" + "Phase Cycle:" + "\t"
            pcyc_table += "["
            dets = self.pcyc[1]
            phases = self.pcyc[0]
            for i in range(0, len(phases)):
                if dets[i] == 1:
                    pcyc_table += "+"
                elif dets[i] == -1:
                    pcyc_table += "-"
                elif dets[i] == 1j:
                    pcyc_table += "+i"
                elif dets[i] == -1j:
                    pcyc_table += "-i"

                if phases[i] == 0:
                    pcyc_table += "(+x)"
                elif phases[i] == np.pi:
                    pcyc_table += "(-x)"
                elif phases[i] == np.pi/2:
                    pcyc_table += "(+y)"
                elif phases[i] == -np.pi/2:
                    pcyc_table += "(-y)"
                elif phases[i] == 3*np.pi/2:
                    pcyc_table += "(-y)"
                pcyc_table += " "
            
            pcyc_table = pcyc_table[:-1]
            pcyc_table += "]\n"
        else:
            pcyc_table = ""

        # Build Footer
        footer = "#" * 79 + "\n" +\
            f"Built by autoDEER Version: {__version__}" + "\n" + "#" * 79

        # Combine All
        string = header + overview_params + param_table + prog_table +\
            pcyc_table + footer

        return string
        

class Parameter:

    def __init__(self, name, value, unit=None, description=None) -> None:
        self.name = name
        self.value = value
        self.unit = unit
        self.description = description
        self.progressive = False
        self.prog = []
        pass

    def add_progression(self, axis_id, axis):
        self.progressive = True
        self.prog.append([axis_id, axis])
    
    def add_step_progression(self, axis_id, length, start, step):
        axis = np.linspace(start, start + step * length, length)
        return self.add_progression(axis_id, axis)


class Detection(Pulse):

    def __init__(self, t, tp) -> None:
        self.t = Parameter("Time", t, "ns", "Start time of pulse")
        self.tp = Parameter("Length", tp, "ns", "Length of the pulse")
        self.Progression = False
        self.pcyc = None
        self.scale = None
        pass


# =============================================================================
# Standard Pulses
# =============================================================================

class RectPulse(Pulse):

    def __init__(self, t, tp, freq, scale) -> None:
        Pulse.__init__(self, t, tp, scale)
        self.freq = Parameter("Freq", freq, "GHz", "Frequency of the Pulse")
        self.Progression = False
        self.pcyc = None
        self._buildFMAM(self.func)
        pass

    def func(self, ax):
        nx = ax.shape[0]
        AM = np.ones(nx)
        FM = np.zeros(nx)
        return AM, FM


# =============================================================================

class HSPulse(Pulse):

    def __init__(self, t, tp, scale, Horder) -> None:
        Pulse.__init__(self, t, tp, scale)
        self.order = Parameter("Order", Horder, None, "Order of the HS Pulse")
        self.Progression = False
        self.pcyc = None

        pass


# =============================================================================

class ChirpPulse(Pulse):

    def __init__(self, t, tp, scale, **kwargs) -> None:
        Pulse.__init__(self, t, tp, scale)

        # Frequency Infomation
        if "BW" in kwargs:
            # Bandwidth + one other
            self.BW = Parameter(
                "Bandwidth", kwargs["BW"], "GHz", "Bandwidth of pulse")
            if "init_freq" in kwargs:
                self.init_freq = Parameter(
                    "f_init", kwargs["init_freq"], "GHz",
                    "Initial frequency of pulse")
            elif "final_freq" in kwargs:
                self.init_freq = Parameter(
                    "f_final", kwargs["final_freq"], "GHz",
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
                self.init_freq.value - self.BW.value, self.final_freq.value,
                nx)
        elif hasattr(self, "init_freq") & hasattr(self, "final_freq"):
            FM = np.linspace(
                self.init_freq.value, self.final_freq.value, nx)

        return AM, FM
        

# =============================================================================

class SincPulse(Pulse):

    def __init__(self, t, tp, scale, freq, order, window=None) -> None:
        Pulse.__init__(self, t, tp, scale)
        self.freq = Parameter("Freq", freq, "GHz", "Frequency of the Pulse")

        self.order = Parameter("Order", order, None, "The sinc pulse order")
        self.window = Parameter(
            "Window", window, None, "The type of window function")

        pass


# =============================================================================
# Composite Pulses
# =============================================================================


class ChorusPulse:

    def __init__(self, t, tp, BW, centre_freq, scale) -> None:
        self.t = t
        self.tp = tp
        self.scale = scale
        pass


# =============================================================================
# Standard Sequences Pulses     !! Move to autoDeer !!
# =============================================================================

def build_HahnEcho(tau, freq):

    pass


# !! Move to autoDeer !!
def build_rectDEER(tau1, tau2, f1, f2, n_pulses=4, nDEER=False):

    pass


def build_AWGDEER(tau1, tau2, f1, f2, n_pulses=4, nDEER=False):

    pass
