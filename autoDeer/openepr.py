import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import hilbert
from scipy.io import savemat
import os
from autodeer import __version__
from autodeer.utils import build_table
import copy
import time
from itertools import product


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

        return fig
        
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

    def add_variable(self, param):
        setattr(self, param.name, param)
        

# =============================================================================


class Interface:

    def __init__(self) -> None:
        pass

    def connect(self) -> None:
        pass

    def acquire_dataset(self) -> dataset:
        pass

    def launch(self, sequence, savename: str):
        pass

    def isrunning(self) -> bool:
        return False

    def terminate(self) -> None:
        pass

    def terminate_at(self, criterion, test_interval=10):
        """Terminates the experiment upon a specific condition being
        satisified. 

        Parameters
        ----------
        criterion : _type_
            The criteria to be tested.
        test_interval : int, optional
            How often should the criteria be tested in minutes, by default 10.
        """

        test_interval_seconds = test_interval * 60
        condition = False

        start_time = time.time()

        while not condition:
            data = self.acquire_dataset()
            try:
                nAvgs = data.nAvgs.value
            except AttributeError:
                print("WARNING: Dataset missing number of averages(nAvgs)!")
                nAvgs = 1
            finally:
                if nAvgs < 1:
                    time.sleep(30)
                    continue

            condition = criterion.test(data)

            if not condition:
                end_time = time.time()

                if (end_time - start_time) < test_interval_seconds:
                    time.sleep(test_interval_seconds - (end_time - start_time))

        self.terminate()
        pass


# =============================================================================
# Super Classes Pulses
# =============================================================================

class Sequence:

    def __init__(
            self, *, name, B, LO, reptime, averages, shots, **kwargs) -> None:
        """Represents an experimental pulse sequence.

        Parameters
        ----------
        name : str
            The name of this pulse sequence
        B : float
            The magnetic field for this sequence in Gauss.
        LO : float
            The central frequency of this sequence. I.e. The frequnecy at which
            a zero offset pulse is at. 
        reptime : float
            The shot repetition time in us.
        averages : int
            The number of scans to be accumulated.
        shots : itn
            The number of shots per point.
        """

        self.pulses = []
        self.num_pulses = len(self.pulses)

        self.B = Parameter(
            "B", B, "Gauss",
            "The static B0 field for the experiment")
        
        self.LO = Parameter(
            "LO", LO, "GHz",
            "The local oscillator frequency.")
        
        self.reptime = Parameter(
            "reptime", reptime, "us",
            "The shot repetition time")
        
        self.averages = Parameter(
            "averages", averages, "None",
            "The number of averages to perform.")
        
        self.shots = Parameter(
            "shots", shots, "None",
            "The number of shots per scan.")

        if "det_window" in kwargs:
            self.det_window = Parameter(
                "det_window", kwargs["det_window"], "None",
                "The length of the default detection gate"
            )
        else:
            self.det_window = Parameter(
                "det_window", 128, "None",
                "The length of the default detection gate"
            )

            
        self.name = name
        pass

    def plot(self) -> None:

        pass

    def addPulse(self, pulse):
        """Adds a pulse to the next position in the sequence.

        Parameters
        ----------
        pulse : Pulse
            The object describing the pulse.
        """
        if type(pulse) == Pulse or issubclass(type(pulse), Pulse):
            self.pulses.append(pulse)
        
        elif type(pulse) == list:
            for el in pulse:
                self.pulses.append(el)
        self.num_pulses = len(self.pulses)
        self._buildPhaseCycle()

    def _estimate_time(self):
        """Calculates the estimated experiment time in seconds.
        """
        acqs = self.averages.value * self.shots.value
        if hasattr(self, 'pcyc_dets'):
            acqs *= self.pcyc_dets.shape[0]
        if hasattr(self, 'progTable'):
            _, pos = np.unique(self.progTable[0], return_index=True)
            for i in pos:
                acqs *= self.progTable[0][i].shape[0]
        time = acqs * self.reptime.value * 1e-6
        return time

    def _buildPhaseCycle(self):
        # Identify pulses which are phase cycled

        pcyc_pulses = []
        pulse_cycles = []
        det_cycles = []
        for ix, pulse in enumerate(self.pulses):
            if pulse.pcyc is not None:
                pcyc_pulses.append(ix)
                # pulse_cycles.append(np.array(pulse.pcyc[0]))
                # det_cycles.append(np.array(pulse.pcyc[1]))
                pulse_cycles.append(pulse.pcyc["Phases"])
                det_cycles.append(pulse.pcyc["DetSigns"])

        self.pcyc_cycles = np.array(list(product(*pulse_cycles)))
        self.pcyc_dets = np.array(list(product(*det_cycles))).prod(axis=1)
        self.pcyc_vars = pcyc_pulses

        # # Build expanded phase cycle
        # func = lambda x: np.arange(0, len(x))
        # map(func, pulse_cycles)
        # n_pulses = len(pulse_cycles)

        # m = list(map(func, pulse_cycles))
        # grids = np.meshgrid(*m, indexing='ij')
        # expanded_cycles = []
        # expanded_dets = []

        # for i in range(0, n_pulses):
        #     expanded_cycles.append(
        #         pulse_cycles[i][grids[i].flatten(order='F')])
        #     expanded_dets.append(det_cycles[i][grids[i].flatten(order='F')])

        # self.pcyc_vars = pcyc_pulses
        # self.pcyc_cycles = np.stack(expanded_cycles)
        # self.pcyc_dets = np.prod(np.stack(expanded_dets), axis=0)
        # self.pcyc_n = self.pcyc_cycles.shape[1]
        return self.pcyc_vars, self.pcyc_cycles, self.pcyc_dets

    def _buildProgTable(self):
        #           parvarid, pulse #, variable, axis
        # progTable = [[], [], [], []]
        progTable = {"EventID": [], "Variable": [], "axis": [],
                     "axID": []}

        for n, pulse in enumerate(self.pulses):
            # Loop over all pulses
            for var_name in vars(pulse):
                var = getattr(pulse, var_name)
                # Loop over all pulse parameters
                if type(var) is Parameter:
                    if var.progressive is True:
                        for i in range(0, len(var.prog)):
                            progTable["axID"].append(var.prog[i][0]) 
                            progTable["EventID"].append(n)
                            progTable["Variable"].append(var_name) 
                            progTable["axis"].append(var.prog[i][1]) 
        
        # Sequence/ Experiment Parameters
        for var_name in vars(self):
            var = getattr(self, var_name)
            if type(var) is Parameter:
                if var.progressive is True:
                    for i in range(0, len(var.prog)):
                        progTable["axID"].append(var.prog[i][0]) 
                        progTable["EventID"].append(None)
                        progTable["Variable"].append(var_name) 
                        progTable["axis"].append(var.prog[i][1]) 

        self.progTable = progTable
        return self.progTable
        
    def addPulseProg(self, pulse_id, variable, axis_id, axis) -> None:
        if pulse_id is None:
            # Assume it is an experiment/sequence parameter
            var = getattr(self, variable)
            var.add_progression(axis_id, axis)
        else:
            var = getattr(self.pulses[pulse_id], variable)
            var.add_progression(axis_id, axis)
        
        self._buildProgTable()
        pass
    
    def addPulsesProg(
            self, pulses, variables, axis_id, axis, multipliers=None) -> None:
        
        if multipliers is None:
            multipliers = np.ones(len(pulses))

        for i, pulse in enumerate(pulses):
            self.addPulseProg(pulse, variables[i], axis_id,
                              axis * multipliers[i])     

        self._buildProgTable()
        pass 

    def isPulseFocused(self):
        test = []
        for pulse in self.pulses:
            test.append(pulse.isPulseFocused())

        if np.array(test).all() == True:
            return True
        else:
            return False

    def isDelayFocused(self):   
        test = []
        for pulse in self.pulses:
            test.append(pulse.isDelayFocused())
        if np.array(test).all() == True:
            return True
        else:
            return False

    def convert(self, *, reduce=True):
        """Converts the current sequence to either pulse focused or delay
        focused depending on the current state

        Parameters
        ----------
        reduce : bool, optional
            Reduce to the smallest number of objects, by default True
        """

        if self.isPulseFocused():
            self._convert_to_delay()
        elif self.isDelayFocused():
            self._convert_to_pulses()

    def _convert_to_delay(self):
        num_pulses = len(self.pulses)

        # create list of pulse timeing
        pulse_times = []
        new_sequence = []

        new_pulse = copy.deepcopy(self.pulses[0])
        new_pulse.t = None
        new_sequence.append(new_pulse)
        pulse_times.append(self.pulses[0].t.value)
        pulse_hash = {0: 0}

        for i in range(1, num_pulses):
            t_cur = self.pulses[i].t.value
            t_prev = self.pulses[i-1].t.value
            static_delay = t_cur-t_prev
            tmp_delay = Delay(tp=static_delay)
            new_sequence.append(tmp_delay)
            new_pulse = copy.deepcopy(self.pulses[i])
            new_pulse.t = None
            new_sequence.append(new_pulse)
            pulse_hash[i] = new_sequence.index(new_pulse)
            pulse_times.append(t_cur)

        new_prog_table = {"EventID": [], "Variable": [], "axis": [],
                          "axID": []}

        def add_prog_table(dic, Element, Variable, axes, ax_id):
            dic["EventID"].append(Element)
            dic["Variable"].append(Variable)
            dic["axis"].append(axes)
            dic["axID"].append(ax_id)
            return dic
        prog_axes = np.unique(self.progTable["axID"])
        for ax_id in prog_axes:
            indexs = np.where(np.array(self.progTable["axID"]) == ax_id)[0]
            for i in indexs:
                pulse_num = self.progTable["EventID"][i]
                pulse_times[pulse_num] = self.progTable["axis"][i]
            for i in range(1, num_pulses):
                diff = np.atleast_1d(pulse_times[i] - pulse_times[i-1])
                new_prog_table = add_prog_table(
                    new_prog_table, int(i*2 - 1), "tp", diff, ax_id)
        
        # Change the pcyc_vars
        new_pcyc_var = []
        for var in self.pcyc_vars:
            new_pcyc_var.append(pulse_hash[var])
        self.pcyc_vars = new_pcyc_var

        self.pulses = new_sequence
        self.progTable = new_prog_table
        return self.pulses

    def _convert_to_pulses(self):
        num_ele = len(self.pulses)
        new_sequence = []

        new_pulse = copy.deepcopy(self.pulses[0])
        new_pulse.t = Parameter("t", 0, "ns", "Start time of pulse")
        new_sequence.append(new_pulse)

        for i in range(1, num_ele):
            pulse = self.pulses[i]
            if type(pulse) is not Delay:
                pos = 0
                for j in range(0, i):
                    pos += self.pulses[j].tp.value
                new_pulse = copy.deepcopy(pulse)
                new_pulse.t = Parameter("t", pos, "ns", "Start time of pulse")
                new_sequence.append(new_pulse)
        self.pulses = new_sequence
        return self.pulses

    def _checkRect(self) -> bool:
        """Checks if all the pulses in the sequence are rectangular.
        """
        test = True

        for pulse in self.pulses:
            if type(pulse) is not RectPulse:
                test = False

        return test

    def __str__(self):

        header = "#" * 79 + "\n" + "AutoDEER Sequence Definition" + \
                 "\n" + "#" * 79 + "\n"

        # Sequence Parameters
        seq_param_string = "Sequence Parameters: \n"
        seq_param_string += "{:<10} {:<12} {:<10} {:<30} \n".format(
            'Name', 'Value', 'Unit', 'Description')

        for param_key in vars(self):
            param = getattr(self, param_key)
            if type(param) is Parameter:
                seq_param_string += "{:<10} {:<12} {:<10} {:<30} \n".format(
                    param.name, param.value, param.unit, param.description)
        
        # Pulses
        pulses_string = "Pulses: \n"

        if self.isPulseFocused():
            # pulses_string += "{:<4} {:<12} {:<8} {:<12} \n".format(
            #     'iD', 't (ns)', 'tp (ns)', 'Type')
            # for i, pulse in enumerate(self.pulses):
            #     pulses_string += "{:<4} {:<12.3E} {:<8} {:<10} \n".format(
            #         i, pulse.t.value, pulse.tp.value, type(pulse).__name__)
            params = ['iD', 't', 'tp', 'scale', 'type']
            params_widths = ["4", "8", "8", "8", "14"]
        elif self.isDelayFocused():
            # pulses_string += "{:<4} {:<8} {:<12} \n".format(
            #     'iD', 'tp (ns)', 'Type')
            # for i, pulse in enumerate(self.pulses):
            #     pulses_string += "{:<4} {:<8} {:<10} \n".format(
            #         i, pulse.tp.value, type(pulse).__name__)
            params = ['iD', 'tp', 'scale', 'type']
            params_widths = ["4", "8", "8", "14"]
        pulses_string += build_table(self.pulses, params, params_widths)
        
        # Progressive elements
        prog_string = "Progression: \n"
        prog_string += "{:<10} {:<10} {:<10} {:<30} \n".format(
            'Pulse', 'Prog. Axis', 'Parameter', 'Step')
        for i in range(0, len(self.progTable["axID"])):
            prog_table = self.progTable
            axis = prog_table["axis"][i]
            if len(np.unique(np.diff(axis))) == 1:
                step = np.unique(np.diff(axis))[0]
            else:
                step = "Var"
            prog_string += "{:<10} {:<10} {:<10} {:<30} \n".format(
                prog_table["EventID"][i], prog_table["axID"][i],
                prog_table["Variable"][i], step)

        footer = "#" * 79 + "\n" +\
            f"Built by autoDEER Version: {__version__}" + "\n" + "#" * 79

        return header + seq_param_string + pulses_string + prog_string + footer

# =============================================================================


class Pulse:

    def __init__(self, *, tp, t, scale=None, flipangle=None, pcyc=None,
                 name=None) -> None:
        """The class for a general pulse.

        Parameters
        ----------
        tp : float
            The pulse length in ns.
        scale : float
            The arbitary experimental pulse amplitude, 0-1.
        t : float, optional
            The pulse start time in ns.

        """
        if t is not None:
            self.t = Parameter("t", t, "ns", "Start time of pulse")
        else:
            self.t = None
        self.tp = Parameter("tp", tp, "ns", "Length of the pulse")
        self.scale = Parameter("scale", scale, None, "Amplitude of pulse")
        self.Progression = False

        self.name = name

        if flipangle is not None:
            self.flipangle = Parameter(
                "flipangle", flipangle, None,
                "The target flip angle of the spins")
        if pcyc is None:
            self.pcyc = [0]
        elif type(pcyc) is dict:
            self._addPhaseCycle(pcyc["phases"], detections=pcyc["dets"])
        else:
            self._addPhaseCycle(pcyc, detections=None)
        pass

    def _addPhaseCycle(self, phases, detections=None):
        if detections is None:
            detections = np.ones(len(phases))
        self.pcyc = {"Phases": list(phases), "DetSigns": list(detections)}
        pass

    def _buildFMAM(self, func):
        self.ax = np.linspace(-self.tp.value/2, self.tp.value/2, 1000)
        self.AM, self.FM = func(self.ax)
        self.complex = self.AM *\
            (np.cos(self.FM*self.ax) + 1j * np.sin(self.FM*self.ax))
        self.FM
        self.AM
        return self.AM, self.FM

    def isDelayFocused(self):
        if self.t is None:
            return True
        else:
            return False
    
    def isPulseFocused(self):
        if self.t is not None:
            return True
        else:
            return False

    def plot(self, pad=1000):
        """Plots the time domain representation of this pulse.

        Parameters
        ----------
        pad : int, optional
            The number of zeros to pad the data with, by default 1000
        """
        dt = self.ax[1] - self.ax[0]
        if self.isPulseFocused():
            tax = self.t.value + self.ax
        else:
            tax = self.ax
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
        header = "#" * 79 + "\n" + "autoDEER Pulse Definition" + \
                 "\n" + "#" * 79 + "\n"
        
        # Build Overviews
        if self.isPulseFocused():
            overview_params = "Time Pos (ns) \tLength (ns) \tScale (0-1)" +\
                "\tProgressive \n"
            
            if type(self) is Detection:
                overview_params += f"{self.t.value}" + "\t\t" + \
                    f"{self.tp.value}" + "\t\t" + "\t\t" +\
                    f"{self.Progression}" + "\n" + "#" * 79 + "\n"        
            else:
                overview_params += f"{self.t.value}" + "\t\t" + \
                    f"{self.tp.value}" + "\t\t" + f"{self.scale.value}" + \
                    "\t\t" + f"{self.Progression}" + "\n" + "#" * 79 + "\n"
        elif self.isDelayFocused():
            overview_params = "Length (ns) \tScale (0-1)" +\
                "\tProgressive \n"
            
            if type(self) is Detection:
                overview_params += f"{self.tp.value}" + "\t\t" + "\t\t" +\
                    f"{self.Progression}" + "\n" + "#" * 79 + "\n"        
            else:
                overview_params += f"{self.tp.value}" + "\t\t" +\
                    f"{self.scale.value}" + "\t\t" + f"{self.Progression}" +\
                    "\n" + "#" * 79 + "\n"

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
            dets = self.pcyc["DetSigns"]
            phases = self.pcyc["Phases"]
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
        """A general parameter.

        Parameters
        ----------
        name : str
            The parameter name
        value : float or int
            The parameter value
        unit : str, optional
            The unit of parameter, by default None. Leave as None if unitless.
        description : str, optional
            A brief description of the parameter, by default None

        Attributes
        ----------
        progressive : bool
            Is the parameter used in any progression or is it constant
        prog : dict
            A dict containing progressive programs for this parameter. This 
            list has two elements. 1) The axis_id"s and 2) the "axis" of 
            values.
        """
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

    def __eq__(self, __o: object) -> bool:
        if type(__o) is not Parameter:
            raise ValueError(
                "Equivalence only works between Parameter classes")
        return self.value == __o.value
        

class Detection(Pulse):

    def __init__(self, *, tp, t) -> None:
        """A general detection pulse.

        Parameters
        ----------
        tp : float
            The **total** time of the detection event. The detection event will
            be symetrical about the centre time. 
        t : float, optional
            The **centre** time of the detection event

        """
        if t is not None:
            self.t = Parameter("Time", t, "ns", "Start time of pulse")
        else:
            self.t = None
        self.tp = Parameter("Length", tp, "ns", "Length of the pulse")
        self.Progression = False
        self.pcyc = None
        self.scale = None
        pass


class Delay(Pulse):

    def __init__(self, *, tp, t=None) -> None:
        
        if t is not None:
            self.t = Parameter("Time", t, "ns", "Start time of pulse")
        else:
            self.t = None
        self.tp = Parameter("Length", tp, "ns", "Length of the pulse")
        self.Progression = False
        self.pcyc = None
        self.scale = None


# =============================================================================
# Standard Pulses
# =============================================================================

class RectPulse(Pulse):

    def __init__(
            self, tp, freq, t=None, flipangle=None, pcyc=None,
            name=None) -> None:
        Pulse.__init__(
            self, tp=tp, t=t, flipangle=flipangle, pcyc=pcyc, name=name)
        self.freq = Parameter("freq", freq, "GHz", "Frequency of the Pulse")
        self.Progression = False
        self._buildFMAM(self.func)
        pass

    def func(self, ax):
        nx = ax.shape[0]
        AM = np.ones(nx)
        FM = np.zeros(nx)
        return AM, FM


# =============================================================================

class HSPulse(Pulse):

    def __init__(self, *, tp, scale, Horder, t) -> None:
        Pulse.__init__(self, tp=tp, scale=scale, t=t)
        self.order = Parameter("Order", Horder, None, "Order of the HS Pulse")
        self.Progression = False
        self.pcyc = None

        pass


# =============================================================================

class ChirpPulse(Pulse):

    def __init__(self, *, tp, scale, t, **kwargs) -> None:
        Pulse.__init__(self, tp=tp, scale=scale, t=t)

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
                self.init_freq = Parameter(
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
                self.init_freq.value - self.BW.value, self.final_freq.value,
                nx)
        elif hasattr(self, "init_freq") & hasattr(self, "final_freq"):
            FM = np.linspace(
                self.init_freq.value, self.final_freq.value, nx)

        return AM, FM
        

# =============================================================================

class SincPulse(Pulse):

    def __init__(self, *, tp, scale, freq, order, t=None, window=None) -> None:
        Pulse.__init__(self, tp=tp, scale=scale, t=t)
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

def build_HahnEcho(tau, pulse_tp, freq, LO, B, scale):
    Hahn_echo = Sequence(
       LO=LO, averages=1, reptime=4e3, shots=200, B=B, name="Hahn_Echo")
    Hahn_echo.addPulse(RectPulse(t=0, tp=pulse_tp, freq=freq, scale=scale))
    Hahn_echo.addPulse(RectPulse(t=tau, tp=pulse_tp*2, freq=freq, scale=scale))
    Hahn_echo.addPulse(Detection(t=tau+pulse_tp*2, tp=512))

    Hahn_echo.addPulsesProg(
        [1, 2],
        ["t", "t"],
        0,
        np.linspace(500, 1000, 21),
        multipliers=[1, 2])
    Hahn_echo._buildProgTable()

    Hahn_echo.pulses[0]._addPhaseCycle([0, np.pi], [1, -1])

    Hahn_echo._buildPhaseCycle()

    return Hahn_echo


def build_FieldSweep(pulse, freq, B, Bsweep=300):

    if type(pulse) is not RectPulse:
        raise RuntimeError("Only rectangular pi/2 pulse are supported here")

    tau = 500

    pulse.pcyc = None
    pulse_pi = copy.deepcopy(pulse)
    pulse_pi.tp.value = pulse.tp.value * 2
    pulse_pi.t.value = tau
    pulse.t.value = 0

    tune = Sequence(
        LO=freq, averages=1, reptime=4e3, shots=100, B=B, name="Tune"
        )

    tune.addPulse(pulse)
    tune.addPulse(pulse_pi)
    tune.addPulse(Detection(t=2*tau + pulse_pi.tp.value, tp=512))

    tune.addPulseProg(
        None,
        'B',
        0,
        np.linspace(B-Bsweep/2, B+Bsweep/2, Bsweep + 1)
    )

    tune._buildProgTable()

    tune.pulses[0]._addPhaseCycle([0, np.pi], [1, -1])

    tune._buildPhaseCycle()

    return tune


# !! Move to autoDeer !!
def build_rectDEER(
        *, tau1, tau2, f1, f2, LO, B, pulse_tp, scale, step=16, n_pulses=4,
        tau3=200):
   
    deadtime = 100

    if n_pulses == 5:
        moving_pulse = [3]
        shift_pulse = 1
        extra_pump = RectPulse(
            t=tau1-tau3, tp=pulse_tp, freq=f2, flipangle=np.pi)
        axis = np.arange(tau1+deadtime, tau2 + 2*tau1 - deadtime, step)

    else:
        moving_pulse = [2]
        extra_pump = None
        shift_pulse = 0
        axis = np.arange(2*tau1-deadtime, tau2 + 2*tau1 - deadtime, step)

    # axis = np.arange(tau1+deadtime, tau2 + tau1 + deadtime, step)
    
    DEER = Sequence(
       LO=LO, averages=1, reptime=4e3, shots=200, B=B, name="autoDEER")
    
    DEER.addPulse(
        RectPulse(t=0, tp=pulse_tp, freq=f1, flipangle=np.pi/2))
    DEER.addPulse(extra_pump)  # pi pump
    DEER.addPulse(RectPulse(t=tau1, tp=pulse_tp, freq=f1, flipangle=np.pi))
    DEER.addPulse(
        RectPulse(
            t=tau1+deadtime, tp=pulse_tp, freq=f2, flipangle=np.pi))  # pi pump
    DEER.addPulse(
        RectPulse(t=2*tau1+tau2, tp=pulse_tp, freq=f1, flipangle=np.pi))  # pi
    DEER.addPulse(Detection(t=2*(tau1+tau2), tp=512))

    DEER.addPulsesProg(
        moving_pulse,
        ["t"],
        0,
        axis)

    DEER._buildProgTable()

    DEER.pulses[1+shift_pulse]._addPhaseCycle(
        [0, np.pi/2, np.pi, -np.pi/2], [1, -1, 1, -1])
    DEER.pulses[2+shift_pulse]._addPhaseCycle(
        [0, np.pi/2, np.pi, -np.pi/2], [1, 1, 1, 1])

    DEER._buildPhaseCycle()
    
    return DEER


def build_AWGDEER(tau1, tau2, f1, f2, n_pulses=4, nDEER=False):

    pass


def tune_pulse(pulse, freq, B):

    if type(pulse) is RectPulse:

        pulse_pi = copy.deepcopy(pulse)
        pulse_pi.tp.value = pulse.tp.value * 2
        pulse_pi.t.value = 500

        tune = Sequence(
            LO=freq, averages=1, reptime=4e3, shots=200, B=B, name="Tune"
        )

        tune.addPulse(pulse)
        tune.addPulse(pulse_pi)
        tune.addPulse(Detection(t=1000 + pulse_pi.tp.value, tp=512))

        tune.addPulsesProg(
            [0, 1],
            ["scale", "scale"],
            0,
            np.linspace(0, 1, 51)
            )
        tune._buildProgTable()

        tune.pulses[0]._addPhaseCycle([0, np.pi], [1, -1])

        tune._buildPhaseCycle()
        
        return tune


def resonator_profile(pulse, freq, gyro):

    if type(pulse) is not RectPulse:
        raise RuntimeError("Only rectangular pi/2 pulse are supported here")
    
    tau0 = 5000
    tau = 500

    hard_pulse = RectPulse(t=0, tp=5.5, freq=0, scale=1.0)
    pulse.pcyc = None
    pulse_pi = copy.deepcopy(pulse)
    pulse_pi.tp.value = pulse.tp.value * 2
    pulse.t.value = tau0
    pulse_pi.t.value = tau0 + tau

    tune = Sequence(
        LO=freq, averages=1, reptime=3e3, shots=50, B=freq/gyro, name="Tune"
        )
    
    tune.addPulse(hard_pulse)
    tune.addPulse(pulse)
    tune.addPulse(pulse_pi)
    tune.addPulse(Detection(t=tau0 + 2*tau + pulse_pi.tp.value, tp=512))

    tune.addPulseProg(
        pulse_id=0,
        variable='tp',
        axis_id=0,
        axis=np.linspace(0, 63, 64)
    )

    tune.addPulsesProg(
        pulses=[0, 1, 2],
        variables=['freq', 'freq', 'freq'],
        axis_id=1,
        axis=np.linspace(-0.4, 0.4, 81)
    )

    tune.addPulseProg(
        pulse_id=None,
        variable='B',
        axis_id=1,
        axis=(np.linspace(-0.4, 0.4, 81)+tune.LO.value)/gyro

    )

    tune._buildProgTable()

    return tune
