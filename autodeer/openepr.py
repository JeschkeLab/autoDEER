import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import hilbert
from scipy.integrate import cumulative_trapezoid
import scipy.fft as fft
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
    Represents an experimental dataset.
    """
    def __init__(self, axes: np.ndarray, data: np.ndarray, params: dict = None,
                 scans: np.ndarray = None) -> None:
        """
        Parameters
        ----------
        axes : list[np.ndarray]
            An array of vectors containing the axes of the dataset
        data : np.ndarray
            The data either as an array of n-dimensions, containing either 
            float or complex values.
        params : dict, optional
            A dictionary of experimental parameters, by default None
        """
        if type(axes) != list:
            self.axes = [self.axes]
        else:
            self.axes = axes
        self.data = data
        self.params = params
        if type(self.axes) is np.ndarray:
            self.dims = self.axes.ndim
        else:
            self.dims = len(self.axes)
        self.scans = scans

        if not np.iscomplexobj(self.data):
            self.data = hilbert(self.data)
        pass

    def plot(self, label=None, **kwargs) -> plt.figure:
        """
        Produces a standard quick graph of the data. This is a line plot for 
        1D data and a heatmap for 2D data.

        Parameters
        ----------
        label : list[str], optional
            A list contains labels. [x_label,z_label], by default None
        
        Optional Kwargs
        ----------
        lines : list[str], optional
            A list of lines to plot. A selection of ["real", "imag","abs"]
            by default all options are plotted.

        Returns
        -------
        plt.figure
            _description_
        """

        if "lines" in kwargs:
            if 'abs' in kwargs["lines"]:
                abs = True
            if 'imag' in kwargs["lines"]:
                imag = True
            if 'abs' in kwargs["lines"]:
                real = True
        else:
            abs=True
            real=True
            imag=True
        if self.dims == 1:
            fig, ax = plt.subplots()
            if abs:
                ax.plot(self.axes[0], np.abs(self.data), label='abs')
            if real:
                ax.plot(self.axes[0], np.real(self.data), label='real')
            if imag:
                ax.plot(self.axes[0], np.imag(self.data), label='imag')
            ax.legend()
            ax.set_ylabel('Signal')
            if label is not None:
                ax.set_xlabel(label[0])
        else:
            raise RuntimeError("Only Single dimension data is supported")

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
    """Represents the interface connection from autoEPR to the spectrometer.
    """

    def __init__(self) -> None:
        pass

    def connect(self) -> None:
        pass

    def acquire_dataset(self) -> dataset:
        """
        Acquires the dataset.
        """
        pass

    def launch(self, sequence, savename: str):
        """Launches the experiment and initialises autosaving.

        Parameters
        ----------
        sequence : Sequence
            The sequence to be launched
        savename : str
            The savename/path for this measurement.
        """
        pass

    def isrunning(self) -> bool:
        return False

    def terminate(self) -> None:
        """
        Terminates the experiment immediately. 
        """
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
    """
    Represents an experimental pulse sequence.
    """

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
        self.progTable = {"EventID": [], "Variable": [], "axis": [],
                     "axID": []}
        pass

    def plot(self) -> None:

        pass

    def plot_pulse_exc(self, FieldSweep=None, ResonatorProfile=None):

        if FieldSweep is not None:
            data = FieldSweep.data
            data /= np.max(np.abs(data))
            ax = FieldSweep.fs_x
            plt.plot(ax, data, label="Field Sweep", color='r')
            
        if ResonatorProfile is not None:
            data = ResonatorProfile.prof_data / ResonatorProfile.prof_data.max()
            plt.plot(ResonatorProfile.prof_frqs-self.LO.value, data, label="Resonator Profile", color='k')


        def check_array_new(list,array, return_pos=False):
            test = True
            for i, ele in enumerate(list):
                if np.allclose(
                    ele,array,rtol=1e-03, atol=1e-03, equal_nan=True):
                    test = False
                    if return_pos:
                        return i
            return test

        axs = []
        fts = []
        labels = []
        for i, pulse in enumerate(self.pulses):
            if type(pulse) == Delay:
                continue
            elif type(pulse) == Detection:
                continue
            ax, ft = pulse._calc_fft(1e4)

            if check_array_new(fts, ft):
                axs.append(ax)
                fts.append(ft)
                labels.append([i])
            else:
                j = check_array_new(fts, ft, return_pos=True)
                labels[j].append(i)

        hatches = ['/', '\\', '|', '-', 'o', '.']
        
        
        for i in range(0,len(axs)):
            ft = np.abs(fts[i])
            ft /= np.max(ft)
            plt.fill(axs[i],ft, 
                     label=f"Pulse id={str(labels[i])[1:-1]}",
                     hatch=hatches.pop(), alpha=0.3)
        
        plt.xlabel("Frequency (GHz)")
        plt.legend()

        # Set correct axis, find mim and max freq
        
        for i, pulse in enumerate(self.pulses):
            if hasattr(pulse, "freq"):
                p_min = pulse.freq.value
                p_max = pulse.freq.value
            elif hasattr(pulse, "init_freq") & hasattr(pulse, "final_freq"):
                p_min = pulse.init_freq.value
                p_max = pulse.final_freq.value
            elif hasattr(pulse, "init_freq") & hasattr(pulse, "BW"):
                p_min = pulse.init_freq.value
                p_max = p_min + pulse.BW.value
            elif hasattr(pulse, "final_freq") & hasattr(pulse, "BW"):
                p_max = pulse.init_freq.value
                p_min = p_max - pulse.BW.value
            # Check p_min < p_max
            if p_min > p_max:
                p_temp = p_min
                p_min = p_max
                p_max = p_temp
            if i == 0:
                min_frq = p_min
                max_frq = p_max
            else:
                if min_frq > p_min:
                    min_frq = p_min
                if max_frq < p_max:
                    max_frq = p_max
        
        min_frq = np.floor(min_frq*100)/100 - 0.2
        max_frq = np.ceil(max_frq*100)/100 + 0.2
        plt.xlim(min_frq, max_frq)
        
        
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
        self._estimate_time()

    def _estimate_time(self):
        """
        Calculates the estimated experiment time in seconds.
        """
        self._buildPhaseCycle()
        acqs = self.averages.value * self.shots.value
        if hasattr(self, 'pcyc_dets'):
            acqs *= self.pcyc_dets.shape[0]
        if hasattr(self, 'progTable'):
            _, pos = np.unique(self.progTable['axID'], return_index=True)
            for i in pos:
                acqs *= self.progTable['axis'][i].shape[0]
        
        time = acqs * self.reptime.value * 1e-6

        self.time = Parameter(name="time", value=f"{(time // 3600):.0f}:{(time % 3600) // 60:.0f}:{(time % 60):.0f}", unit="HH:MM:SS",
        description="Estimated sequence run time")
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
        """Adds a single progessive pulse element to the sequence.

        It is strongly recomeneded that the axis is generated using `np.arange`
        over `np.linspace`. Most spectrometers do not like inputs with a high
        number of decimal places. `np.arange` allows control of the step size
        reducing any risk of conflict with the interface.

        Parameters
        ----------
        pulse_id : int
            The iD of the moving pulse element. If it is a sequence parameter, 
            such as `B` field, then `None` should be given.
        variable : str
            The name of the parameter.
        axis_id : int
            The iD of the axis. 
        axis : np.ndarray
            An array containing how the pulse element changes.
        """
        if pulse_id is None:
            # Assume it is an experiment/sequence parameter
            var = getattr(self, variable)
            var.add_progression(axis_id, axis)
        else:
            var = getattr(self.pulses[pulse_id], variable)
            var.add_progression(axis_id, axis)
        
        self._buildProgTable()
        self._estimate_time()

        pass
    
    def addPulsesProg(
            self, pulses, variables, axis_id, axis, multipliers=None) -> None:
        """Adds a multiple progessive pulse element to the sequence. This is
        very useful when multiple elements are changing at the same time with
        a constant relationship.

        It is strongly recomeneded that the axis is generated using `np.arange`
        over `np.linspace`. Most spectrometers do not like inputs with a high
        number of decimal places. `np.arange` allows control of the step size
        reducing any risk of conflict with the interface.


        Parameters
        ----------
        pulses : list[int]
            A list of pulse iDs that are changing
        variables : list[str]
            A list of variables that are changing
        axis_id : int
            The iD of the axis 
        axis : np.ndaaray
            An array containing how the pulse element changes
        multipliers : list or np.ndaaray, optional
            How the different variable are proportional to each other, by
            default None. If `None` is specified then it is assumed that there
            is a 1:1 relationship. The axs gives the values for the first
            element.
        """
        
        if multipliers is None:
            multipliers = np.ones(len(pulses))

        for i, pulse in enumerate(pulses):
            self.addPulseProg(pulse, variables[i], axis_id,
                              axis * multipliers[i])     

        self._buildProgTable()
        pass 

    def isPulseFocused(self):
        """
        Is the sequence expressed to contain only pulses and no delays?
        """
        test = []
        for pulse in self.pulses:
            test.append(pulse.isPulseFocused())

        if np.array(test).all() == True:
            return True
        else:
            return False

    def isDelayFocused(self):  
        """
        Is the sequence expressed to contain both pulses and delays?
        """ 
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
                if diff.shape[0] > 1:
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

        header = "#" * 79 + "\n" + "autoDEER Sequence Definition" + \
                 "\n" + "#" * 79 + "\n"

        # Sequence Parameters
        seq_param_string = "Sequence Parameters: \n"
        seq_param_string += "{:<10} {:<12} {:<10} {:<30} \n".format(
            'Name', 'Value', 'Unit', 'Description')

        for param_key in vars(self):
            param = getattr(self, param_key)
            if type(param) is Parameter:
                if type(param.value) is str:
                    seq_param_string += "{:<10} {:<12} {:<10} {:<30} \n".format(
                        param.name, param.value, param.unit, param.description)
                else:
                    seq_param_string += "{:<10} {:<12.5g} {:<10} {:<30} \n".format(
                        param.name, param.value, param.unit, param.description)
        
        # Pulses
        pulses_string = "\nEvents (Pulses, Delays, etc...): \n"

        if self.isPulseFocused():
            # pulses_string += "{:<4} {:<12} {:<8} {:<12} \n".format(
            #     'iD', 't (ns)', 'tp (ns)', 'Type')
            # for i, pulse in enumerate(self.pulses):
            #     pulses_string += "{:<4} {:<12.3E} {:<8} {:<10} \n".format(
            #         i, pulse.t.value, pulse.tp.value, type(pulse).__name__)
            params = ['iD', 't', 'tp', 'scale', 'type', 'Phase Cycle']
            params_widths = ["4", "8", "8", "8", "14", "40"]
        elif self.isDelayFocused():
            # pulses_string += "{:<4} {:<8} {:<12} \n".format(
            #     'iD', 'tp (ns)', 'Type')
            # for i, pulse in enumerate(self.pulses):
            #     pulses_string += "{:<4} {:<8} {:<10} \n".format(
            #         i, pulse.tp.value, type(pulse).__name__)
            params = ['iD', 'tp', 'scale', 'type']
            params_widths = ["4", "8", "8", "14"]
        pulses_string += build_table(self.pulses, params, params_widths)

        def print_event_id(i):
            if prog_table["EventID"][i] is None:
                return "Seq"
            else:
                return str(prog_table["EventID"][i])

        def get_unit(i):
            pulse_num = prog_table["EventID"][i]
            if pulse_num is None:
                param = getattr(self, prog_table["Variable"][i])
            else:
                param = getattr(
                    self.pulses[pulse_num], prog_table["Variable"][i])
            
            if param.unit is None:
                return "None"
            else:
                return param.unit


        def test_unique_step(array):
            diffs = np.diff(array)-np.diff(array)[0]
            return np.isclose(diffs,np.zeros(diffs.shape)).all()

        # Progressive elements
        prog_string = "\nProgression: \n"
        if len(self.progTable["axID"]) >= 1:
            prog_string += "{:<10} {:<10} {:<10} {:<10} {:<10} {:<10} \n".format(
                'Pulse', 'Prog. Axis', 'Parameter', 'Step', 'Dim', 'Unit')
        for i in range(0, len(self.progTable["axID"])):
            prog_table = self.progTable
            axis = prog_table["axis"][i]
            if test_unique_step(axis):
                step = np.unique(np.diff(axis))[0]
                fstring = "{:<10} {:<10} {:<10} {:<10.5g} {:<10} {:<10} \n"
            else:
                step = "Var"
                fstring = "{:<10} {:<10} {:<10} {:<10} {:<10} {:<10} \n"
             
            prog_string += fstring.format(
                print_event_id(i), prog_table["axID"][i],
                prog_table["Variable"][i], step,
                prog_table["axis"][i].shape[0], get_unit(i))

        footer = "#" * 79 + "\n" +\
            f"Built by autoDEER Version: {__version__}" + "\n" + "#" * 79

        return header + seq_param_string + pulses_string + prog_string + footer

    def copy(self):
        return copy.deepcopy(self)

# =============================================================================


class Pulse:
    """
    Represents a general experimental pulse.
    """

    def __init__(self, *, tp, t, scale=None, flipangle=None, pcyc=[0],
                 name=None, **kwargs) -> None:
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
            self.pcyc = None
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
        self.ax = np.linspace(-self.tp.value/2, self.tp.value/2, 1001)
        dt = self.ax[1]-self.ax[0]
        self.AM, self.FM = func(self.ax)
        FM_arg = 2*np.pi*cumulative_trapezoid(self.FM, initial=0) * dt
        self.complex = self.AM * (np.cos(FM_arg) +1j* np.sin(FM_arg))
        self.FM
        self.AM
        return self.AM, self.FM

    def isDelayFocused(self):
        """
        Does the pulse contain a specified time, `t`?

        If so then it is not delay focused.
        """
        if self.t is None:
            return True
        else:
            return False
    
    def isPulseFocused(self):
        """
        Does the pulse contain a specified time, `t`?

        If so then it is delay focused.
        """
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
    
    def _calc_fft(self, pad=10000):
        self._buildFMAM(self.func)
        data = np.pad(self.complex, int(pad))
        pulse_fft = fft.fftshift(fft.fft(data))
        n = data.shape[0]
        dt = self.ax[1] -self.ax[0]
        axis_fft = fft.fftshift(fft.fftfreq(n, dt))

        return axis_fft, pulse_fft
    
    def plot_fft(self):
                
        ax, ft = self._calc_fft(1e4)
        plt.plot(ax, np.real(ft), label="Re")
        plt.plot(ax, np.imag(ft), label="Im")
        plt.plot(ax, np.abs(ft), label="Im")
        plt.legend()
        plt.xlabel("Frequency (GHz)")
        plt.ylabel("Amplitude (A.U.)")
        pass

    def _pcyc_str(self):
    
        if self.pcyc is not None:
            pcyc_table = "["
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
            pcyc_table += "]"
        else:
            pcyc_table = ""

        return pcyc_table

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
                if attr.name == "flipangle":
                    if attr.value == np.pi:
                        value = "π"
                    elif attr.value == np.pi/2:
                        value = "π/2"
                    else:
                        value = f"{attr.value:>5.5g}"
                else:
                    value = f"{attr.value:>5.5g}"

                param_table += f"{attr.name} \t\t {value} \t" +\
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
            pcyc_table += self._pcyc_str()
            pcyc_table += "\n"
        else:
            pcyc_table = ""

        # Build Footer
        footer = "#" * 79 + "\n" +\
            f"Built by autoDEER Version: {__version__}" + "\n" + "#" * 79

        # Combine All
        string = header + overview_params + param_table + prog_table +\
            pcyc_table + footer

        return string

    def copy(self, **kwargs):
        """Creates a deep-copy of the pulse. I.e. Every parameter object is
        re-created at another memory space.

        Parameter can be chaged at this stage by adding them as keyword-
        arguments (kwargs).

        Returns
        -------
        Pulse
            A deep copy of the pulse
        """

        new_pulse = copy.deepcopy(self)

        for arg in kwargs:
            if (hasattr(new_pulse,arg)) and (getattr(new_pulse,arg) is not None):
                attr = getattr(new_pulse,arg)
                if type(attr) is Parameter:
                    attr.value = kwargs[arg]                    
                else:
                    attr = kwargs[arg]
            elif arg == "t":
                setattr(new_pulse, arg,
                        Parameter("t", kwargs[arg], "ns",
                                  "Start time of pulse"))
            elif arg == "scale":
                setattr(new_pulse, arg,
                        Parameter("scale", kwargs[arg], None,
                                  "Amplitude of pulse"))
            elif arg == "flipangle":
                setattr(new_pulse,arg,
                        Parameter("flipangle", kwargs[arg], None,
                                  "The target flip angle of the spins"))
            else:
                setattr(new_pulse,arg,kwargs[arg])
        
        return new_pulse
        

class Parameter:
    """
    Represents a sequence or pulse parameter.
    """

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
    
    def remove_progression(self):
        """Removes progression from parameter object
        """
        self.progressive = False
        self.prog = []
        pass

    def __eq__(self, __o: object) -> bool:
        if type(__o) is not Parameter:
            raise ValueError(
                "Equivalence only works between Parameter classes")
        return self.value == __o.value
    
    def copy(self):
        return copy.deepcopy(self)
        

class Detection(Pulse):
    """
    Represents a detection pulse.
    """

    def __init__(self, *, tp, t=None, freq=0,**kwargs) -> None:
        """A general detection pulse.

        Parameters
        ----------
        tp : float
            The **total** time of the detection event. The detection event will
            be symetrical about the centre time. 
        t : float, optional
            The **centre** time of the detection event
        freq: float, optional
            The detection frequency, not all spectrometer support this 
            functionality, by default 0 MHz

        """
        if t is not None:
            self.t = Parameter("t", t, "ns", "Start time of pulse")
        else:
            self.t = None
        self.tp = Parameter("tp", tp, "ns", "Length of the pulse")
        self.freq = Parameter(
            "freq", freq, "MHz", "The detection frequency offset from the LO")
        self.Progression = False
        self.pcyc = None
        self.scale = None
        pass


class Delay(Pulse):
    """
    Represents a inter-pulse delay pulse.
    """

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
    """
    Represents a CHORUS composite pulse.
    """
    def __init__(self, t, tp, BW, centre_freq, scale) -> None:
        self.t = t
        self.tp = tp
        self.scale = scale
        pass

