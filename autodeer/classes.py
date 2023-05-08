import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import hilbert
from scipy.integrate import cumulative_trapezoid
import scipy.fft as fft
from scipy.io import savemat
import os
from autodeer import __version__
from autodeer.utils import build_table, sop
import copy
import time
from itertools import product
import numbers
import uuid


# =============================================================================
class Dataset:
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
            self.axes = [axes]
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

    def acquire_dataset(self) -> Dataset:
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

    def terminate_at(self, criterion, test_interval=10, verbosity=0):
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
        last_scan = 0


        while not condition:
            start_time = time.time()
            data = self.acquire_dataset()
            try:
                nAvgs = data.nAvgs.value
            except AttributeError:
                print("WARNING: Dataset missing number of averages(nAvgs)!")
                nAvgs = 1
            finally:
                if nAvgs < 1:
                    time.sleep(30)  # Replace with single scan time
                    continue
                elif nAvgs <= last_scan:
                    time.sleep(30)
                    continue    
            last_scan = nAvgs
            if verbosity > 0:
                print("Testing")
            condition = criterion.test(data, verbosity)

            if not condition:
                if not self.isrunning():
                    msg = "Experiments has finished before criteria met."
                    raise RuntimeError(msg)
                end_time = time.time()
                if (end_time - start_time) < test_interval_seconds:
                    if verbosity > 0:
                        print("Sleeping")
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


        if isinstance(B, Parameter):
            self.B = B.copy()
        else:
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
                     "axID": [], "uuid": [], "reduce": []}
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

    def evolution(self, params, reduce=[]):
        axes_uuid = [param.uuid for param in params]
        reduce_uuid = [param.uuid for param in reduce]


        # Build the ProgTable
        progTable = {"EventID": [], "Variable": [], "axis": [],
                     "axID": [], "uuid": [], "reduce": []}
        
        for n, pulse in enumerate(self.pulses):
            table = pulse.build_table()
            for i in range(len(table["uuid"])):
                if table["uuid"][i] in axes_uuid:
                    progTable["axID"].append(axes_uuid.index(table["uuid"][i]))
                    progTable["uuid"].append(table["uuid"][i]) 
                    progTable["EventID"].append(n)
                    progTable["Variable"].append(table["Variable"][i])
                    progTable["axis"].append(table["axis"][i])
                    if table["uuid"][i] in reduce_uuid:
                        progTable["reduce"].append(True)
                    else:
                        progTable["reduce"].append(False)

        for var_name in vars(self):
            var = getattr(self, var_name)
            if type(var) is Parameter:
                if not var.is_static():
                    for i in range(len(var.axis)):
                        if var.axis[i]["uuid"] in axes_uuid:
                            progTable["axID"].append(axes_uuid.index(var.axis[i]["uuid"]))
                            progTable["EventID"].append(None)
                            progTable["Variable"].append(var_name) 
                            progTable["axis" ].append(var.axis[i]["axis"])
                            progTable["uuid"].append(var.axis[i]["uuid"]) 
                            if var.axis[i]["uuid"] in reduce_uuid:
                                progTable["reduce"].append(True)
                            else:
                                progTable["reduce"].append(False)
        self.progTable = progTable
        return self.progTable
        
    def _buildProgTable(self):
        #           parvarid, pulse #, variable, axis
        # progTable = [[], [], [], []]
        progTable = {"EventID": [], "Variable": [], "axis": [],
                     "axID": []}

        for n, pulse in enumerate(self.pulses):
            # Loop over all pulses
            # for var_name in vars(pulse):
            #     var = getattr(pulse, var_name)
            #     # Loop over all pulse parameters
            #     if type(var) is Parameter:
            #         if var.progressive is True:
            #             for i in range(0, len(var.prog)):
            #                 progTable["axID"].append(var.prog[i][0]) 
            #                 progTable["EventID"].append(n)
            #                 progTable["Variable"].append(var_name) 
            #                 progTable["axis"].append(var.prog[i][1]) 

            table = pulse.build_table()
            n_table = len(table["uuid"])
            progTable["axID"] += table["uuid"]
            progTable["EventID"] += [n] * n_table
            progTable["Variable"] += table["Variable"] 
            progTable["axis"] += table["axis"]  

        # Sequence/ Experiment Parameters
        for var_name in vars(self):
            var = getattr(self, var_name)
            if type(var) is Parameter:
                if not var.is_static():
                    for i in range(0, len(var.prog)):
                        # progTable["axID"].append(var.prog[i][0]) 
                        # progTable["EventID"].append(None)
                        # progTable["Variable"].append(var_name) 
                        # progTable["axis"].append(var.prog[i][1]) 
                        progTable["axID"] += var.ax_id
                        progTable["EventID"].append(None)
                        progTable["Variable"].append(var_name) 
                        progTable["axis" ]+= var.axis
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

    def __init__(self, *, tp, t=None, scale=None, flipangle=None, pcyc=[0],
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
            if isinstance(t,numbers.Number):
                self.t = Parameter("t", t, "ns", "Start time of pulse")
            elif isinstance(t,Parameter):
                self.t = t
        else:
            self.t = None

        if isinstance(tp,numbers.Number):
            self.tp = Parameter("tp", tp, "ns", "Length of the pulse")
        elif isinstance(tp,Parameter):
            self.tp = tp

        if isinstance(scale,numbers.Number):
            self.scale = Parameter("scale", scale, None, "Amplitude of pulse")
        elif isinstance(scale,Parameter):
            self.scale = scale


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
    

    def build_table(self):
        progTable = {"Variable": [], "axis": [],
                "uuid": []}
        for arg in vars(self):
            attr = getattr(self,arg)
            if type(attr) is Parameter and not attr.is_static():
                for i, axes in enumerate(attr.axis):
                    progTable["Variable"].append(arg)
                    progTable["axis"].append(axes["axis"])
                    progTable["uuid"].append(axes["uuid"])

        return progTable

    
    def is_static(self):

        for arg in vars(self):
            attr = getattr(self,arg)
            if type(attr) is Parameter and not attr.is_static():
                return False
        return True

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

    def exciteprofile(self, freqs=None):
        """Excitation profile

        Generates the exciatation profiles for this pulse.

        This function is ported from EasySpin 
        (https://easyspin.org/easyspin/documentation/sop.html) [1-2],
        and based upon the method from Gunnar Jeschke, Stefan Pribitzer and
        Andrin Doll[3].

        References:
        +++++++++++
        [1] Stefan Stoll, Arthur Schweiger
        EasySpin, a comprehensive software package for spectral simulation and analysis in EPR
        J. Magn. Reson. 178(1), 42-55 (2006)
        
        [2] Stefan Stoll, R. David Britt
        General and efficient simulation of pulse EPR spectra
        Phys. Chem. Chem. Phys. 11, 6614-6625 (2009)

        [3] Jeschke, G., Pribitzer, S. & DollA. 
        Coherence Transfer by Passage Pulses in Electron Paramagnetic 
        Resonance Spectroscopy. 
        J. Phys. Chem. B 119, 13570-13582 (2015)


        Parameters
        ----------
        
        freqs: np.ndarray, optional
            The frequency axis. Caution: A larger number of points will linearly
            increase computation time.
        
        Returns
        -------

        Mx: np.ndarray
            The magentisation in the X direction.
        My: np.ndarray
            The magentisation in the Y direction.
        Mz: np.ndarray
            The magentisation in the Z direction.
        """
        eye = np.identity(2)
        # offsets = np.arange(-2*BW,2*BW,0.002)

        if freqs is None:
            fxs, fs = self._calc_fft()
            fs = np.abs(fs)
            points = np.argwhere(fs>(fs.max()*0.5))[:,0]
            BW = fxs[points[-1]] - fxs[points[0]]
            c_freq = np.mean([fxs[points[-1]], fxs[points[0]]])
            offsets = np.linspace(-2*BW, 2*BW, 100) + c_freq
        else:
            offsets = freqs

        t = self.ax
        nOffsets = offsets.shape[0]
        amp_factor = self.flipangle.value / (2 * np.pi * np.trapz(self.AM,t))
        ISignal = np.real(self.complex) * amp_factor
        QSignal = np.imag(self.complex) * amp_factor
        npoints = ISignal.shape[0]
        Sx = sop([1/2],"x")
        Sy = sop([1/2],"y")
        Sz = sop([1/2],"z")
        Density0 = -Sz
        Mag = np.zeros((3,nOffsets), dtype=np.complex128)
        for iOffset in range(0,nOffsets):
            UPulse = eye
            Ham0 = offsets[iOffset]*Sz

            for it in range(0,npoints):
                Ham = ISignal[it] * Sx + QSignal[it] * Sy + Ham0
                M = -2j * np.pi * (t[1]-t[0]) * Ham
                q = np.sqrt(M[0,0] ** 2 - np.abs(M[0,1]) ** 2)
                if np.abs(q) < 1e-10:
                    dU = eye + M
                else:
                    dU = np.cosh(q)*eye + (np.sinh(q)/q)*M
                UPulse = dU @ UPulse

            Density = UPulse @ Density0 @ UPulse.conjugate().T
            Mag[0, iOffset] = -2 * (Sx @ Density.T).sum().real
            Mag[1, iOffset] = -2 * (Sy @ Density.T).sum().real
            Mag[2, iOffset] = -2 * (Sz @ Density.T).sum().real
        
        return Mag[0,:], Mag[1,:], Mag[2,:] 


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
                "\t Static \n"
            
            if type(self) is Detection:
                overview_params += f"{self.t.value}" + "\t\t" + \
                    f"{self.tp.value}" + "\t\t" + "\t\t" +\
                    f"{self.is_static()}" + "\n" + "#" * 79 + "\n"        
            else:
                overview_params += f"{self.t.value}" + "\t\t" + \
                    f"{self.tp.value}" + "\t\t" + f"{self.scale.value}" + \
                    "\t\t" + f"{self.is_static()}" + "\n" + "#" * 79 + "\n"
        elif self.isDelayFocused():
            overview_params = "Length (ns) \tScale (0-1)" +\
                "\t Static \n"
            
            if type(self) is Detection:
                overview_params += f"{self.tp.value}" + "\t\t" + "\t\t" +\
                    f"{self.is_static()}" + "\n" + "#" * 79 + "\n"        
            else:
                overview_params += f"{self.tp.value}" + "\t\t" +\
                    f"{self.scale.value}" + "\t\t" + f"{self.is_static()}" +\
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
                elif attr.value is None:
                    value = "N/A"
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

    def copy(self, clear=False, **kwargs):
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

        if clear:
            for arg in vars(new_pulse):
                attr = getattr(new_pulse,arg)
                if type(attr) is Parameter and not getattr(new_pulse,arg).is_static():
                    getattr(new_pulse,arg).remove_dynamic()


        for arg in kwargs:
            if (hasattr(new_pulse,arg)) and (getattr(new_pulse,arg) is not None):
                attr = getattr(new_pulse,arg)
                if type(attr) is Parameter:
                    attr.value = kwargs[arg]     
                elif arg == "pcyc":
                    new_pcyc = kwargs[arg]
                    if attr is None:
                        new_pulse.pcyc = None
                    elif type(new_pcyc) is dict:
                        new_pulse._addPhaseCycle(new_pcyc["phases"], detections=new_pcyc["dets"])
                    else:
                        new_pulse._addPhaseCycle(new_pcyc, detections=None)               
                else:
                    attr = kwargs[arg]
            elif arg == "t":
                if type(kwargs[arg]) is Parameter:
                    kwargs[arg].name = "t"
                    setattr(new_pulse, arg, kwargs[arg])
                elif isinstance(kwargs[arg], numbers.Number):
                    setattr(new_pulse, arg,
                            Parameter("t", kwargs[arg], "ns",
                                    "Start time of pulse"))
                else:
                    raise ValueError("new parameters must be either numeric or Parameters")
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

    def __init__(self, name, value, unit=None, description=None, **kwargs) -> None:
        """A general parameter.

        Parameters
        ----------
        name : str
            The parameter name
        value : float or int
            The parameter value, eithe initial or static
        unit : str, optional
            The unit of parameter, by default None. Leave as None if unitless.
        description : str, optional
            A brief description of the parameter, by default None    
        axis : np.ndarray, optional
            The difference from the intial value for each position in a 
            dynamic axis. Can be n-dimensional, by default None.
        ax_id : list, optional 
             
    

        Attributes
        ----------
        progressive : bool
            Is the parameter used in any progression or is it constant
        prog : dict
            A dict containing progressive programs for this parameter. This 
            list has two elements. 1) The axis_id"s and 2) the "axis" of 
            values.

        Parameter Arthimatic
        --------------------



        Examples
        --------
        Creating a static parameter
        ```
        Par1 = Parameter(
            name="Par1", value=10, unit="us", description="The first parameter")
        ```
        Creating a dynamic parameter
        ```
        ```
        Making a dynamic parameter
        ```
        ```

        Adding a parameter and a number:
        ```
        Par1 = Parameter(
            name="Par1", value=10, unit="us", description="The first parameter")
        Par2 = Par1 + 2
        """

        self.name = name
        self.value = value
        self.unit = unit
        self.description = description
        self.axis = []
        self.ax_id = []
        if "link" in kwargs:
            if not isinstance(kwargs["link"], Parameter):
                raise ValueError("The linked parameter must be a Parmater object")
            self.uuid = kwargs["link"].uuid
        else:
            self.uuid = uuid.uuid1()

        if "step" in kwargs:
            step = kwargs["step"]
            dim = kwargs["dim"]
            if "axis_id" in kwargs:
                axis_id = kwargs["axis_id"]
            else:
                axis_id = 0
            if "start" in kwargs:
                start = kwargs["start"]
            else:
                start = 0
            axis = np.arange(start=0, stop= dim*step+start,step=step)
            self.add_axis(axis=axis,axis_id=axis_id)
        


        pass

    def add_axis(self, axis_id, axis):
        # if self.axis == []:
        #     self.axis.append(np.array(axis))
        #     self.ax_id.append(axis_id)
        self.axis.append({"axis":axis, "uuid":self.uuid})

    def remove_dynamic(self):
        self.axis = []
        self.ax_id =[]     
    
    def is_static(self) -> bool:
        if self.axis == []:
            return True
        # elif self.ax_id == []:
        #     return True
        else:
            return False

    def __eq__(self, __o: object) -> bool:
        if type(__o) is not Parameter:
            raise ValueError(
                "Equivalence only works between Parameter classes")
        return self.value == __o.value
    
    def __add__(self, __o:object):

        if type(__o) is Parameter:
            if self.unit != __o.unit:
                raise RuntimeError("Both parameters must have the same unit")
            new_value = self.value + __o.value
            new_name = f"{self.name} + {__o.name}"
            new_description  = new_name
            new_parameter = Parameter(
                name=new_name, value=new_value, unit=self.unit,
                description=new_description)
            if not self.is_static():
                if not __o.is_static():
                    # Dynamic parmaters can only be summed and multiplied if the axis has the same uuid. I.e. they were linked when created or are deriratives of each other. 
                    new_ax_id = []
                    new_axis = []
                    # a_ax_ids:list = self.ax_id
                    a_ax_ids:list = [self.axis[i]["uuid"] for i in range(len(self.axis))]
                    # b_ax_ids:list = __o.ax_id
                    b_ax_ids:list = [__o.axis[i]["uuid"] for i in range(len(__o.axis))]
                    ab_ax_ids = list(set(a_ax_ids + b_ax_ids))
                    for id in ab_ax_ids:
                        if id not in b_ax_ids: # I.e. only in A
                            a_index = a_ax_ids.index(id)
                            new_axis.append(self.axis[a_index])
                            new_ax_id.append(id)
                        elif id not in a_ax_ids: # I.e. only in B
                            b_index = b_ax_ids.index(id)
                            new_axis.append(__o.axis[b_index])
                            new_ax_id.append(id)
                        else: # in both
                            a_index = a_ax_ids.index(id)
                            b_index = b_ax_ids.index(id)
                            b_ax_ids.remove(id)
                            new_axis.append({"axis": self.axis[a_index]["axis"] + __o.axis[b_index]["axis"], "uuid": self.uuid})
                            new_ax_id.append(id)
                else:
                    new_axis = self.axis
                    new_ax_id = self.ax_id

            else:
                if not __o.is_static():
                    new_axis = __o.axis
                    new_ax_id = __o.ax_id
                else:
                    new_axis = []
                    new_ax_id = []

            new_parameter.axis = new_axis
            new_parameter.ax_id = new_ax_id

            return new_parameter

        elif isinstance(__o, numbers.Number):
            new_value = self.value + __o
            new_name = f"{self.name} + {__o}"
            new_parameter = Parameter(
                name=new_name, value=new_value, unit=self.unit)
            if not self.is_static():
                new_axis = self.axis
                new_ax_id = self.ax_id
                new_parameter.axis = new_axis
                new_parameter.ax_id = new_ax_id
            return new_parameter
        
        elif isinstance(__o, np.ndarray):
            if self.axis.shape != __o.shape:
                raise RuntimeError(
                    "Both parameters axis and the array must have the same shape")
    
    def __sub__(self, __o:object):
        
        if type(__o) is Parameter:
            if self.unit != __o.unit:
                raise RuntimeError("Both parameters must have the same unit")
            new_value = self.value - __o.value
            new_name = f"{self.name} - {__o.name}"
            new_parameter = Parameter(
                name=new_name, value=new_value, unit=self.unit)
            if not self.is_static():
                if not __o.is_static():
                    # Dynamic parmaters can only be summed and multiplied if the axis has the same uuid. I.e. they were linked when created or are deriratives of each other. 
                    new_ax_id = []
                    new_axis = []
                    # a_ax_ids:list = self.ax_id
                    a_ax_ids:list = [self.axis[i]["uuid"] for i in range(len(self.axis))]
                    # b_ax_ids:list = __o.ax_id
                    b_ax_ids:list = [__o.axis[i]["uuid"] for i in range(len(self.__o))]
                    ab_ax_ids = list(set(a_ax_ids + b_ax_ids))
                    for id in ab_ax_ids:
                        if id not in b_ax_ids: # I.e. only in A
                            a_index = a_ax_ids.index(id)
                            new_axis.append({"axis": self.axis[a_index], "uuid": self.uuid})
                            new_ax_id.append(id)
                        elif id not in a_ax_ids: # I.e. only in B
                            b_index = b_ax_ids.index(id)
                            new_axis.append({"axis": __o.axis[b_index], "uuid": __o.uuid})
                            new_ax_id.append(id)
                        else: # in both
                            a_index = a_ax_ids.index(id)
                            b_index = b_ax_ids.index(id)
                            b_ax_ids.remove(id)
                            new_axis.append({"axis": self.axis[a_index] - __o.axis[b_index], "uuid": self.uuid})
                            new_ax_id.append(id)
                else:
                    new_axis = self.axis
                    new_ax_id = self.ax_id

            else:
                if not __o.is_static():
                    new_axis = __o.axis
                    new_ax_id = __o.ax_id
                else:
                    new_axis = []
                    new_ax_id = []

            new_parameter.axis = new_axis
            new_parameter.ax_id = new_ax_id

            return new_parameter

        elif isinstance(__o, numbers.Number):
            new_value = self.value - __o
            new_name = f"{self.name} - {__o}"
            new_parameter = Parameter(
                name=new_name, value=new_value, unit=self.unit)
            if self.axis is not None:
                new_axis = self.axis 
                new_ax_id = self.ax_id
                new_parameter.axis = new_axis
                new_parameter.ax_id = new_ax_id
            return new_parameter
        
        elif isinstance(__o, np.ndarray):
            if self.axis.shape != __o.shape:
                raise RuntimeError(
                    "Both parameters axis and the array must have the same shape")

    def __mul__(self, __o:object):
        if type(__o) is Parameter:
            if self.unit != __o.unit:
                raise RuntimeError("Both parameters must have the same unit")
            # if not __o.is_static():
            #     raise RuntimeError("Multiplictaion of two dynamic parameters is not supported")
            new_value = self.value * __o.value
            new_name = f"{self.name} * {__o.name}"
            new_parameter = Parameter(
                name=new_name, value=new_value, unit=self.unit)
            # if self.axis is not None:
            #     new_axis =  [np.array([item * __o.value for item in axis]) for axis in self.axis ]
            #     new_ax_id = self.ax_id
            #     new_parameter.axis = new_axis
            #     new_parameter.ax_id = new_ax_id
            # return new_parameter
            if not self.is_static():
                if not __o.is_static():
                    # Dynamic parmaters can only be summed and multiplied if the axis has the same uuid. I.e. they were linked when created or are deriratives of each other. 
                    new_ax_id = []
                    new_axis = []
                    # a_ax_ids:list = self.ax_id
                    a_ax_ids:list = [self.axis[i]["uuid"] for i in range(len(self.axis))]
                    # b_ax_ids:list = __o.ax_id
                    b_ax_ids:list = [__o.axis[i]["uuid"] for i in range(len(self.__o))]
                    ab_ax_ids = list(set(a_ax_ids + b_ax_ids))
                    for id in ab_ax_ids:
                        if id not in b_ax_ids: # I.e. only in A
                            a_index = a_ax_ids.index(id)
                            new_axis.append({"axis": self.axis[a_index], "uuid": self.uuid})
                            new_ax_id.append(id)
                        elif id not in a_ax_ids: # I.e. only in B
                            b_index = b_ax_ids.index(id)
                            new_axis.append({"axis": __o.axis[b_index], "uuid": __o.uuid})
                            new_ax_id.append(id)
                        else: # in both
                            a_index = a_ax_ids.index(id)
                            b_index = b_ax_ids.index(id)
                            b_ax_ids.remove(id)
                            new_axis.append({"axis": self.axis[a_index] * __o.axis[b_index], "uuid": self.uuid})
                            new_ax_id.append(id)
                else:
                    new_axis = self.axis
                    new_ax_id = self.ax_id

            else:
                if not __o.is_static():
                    new_axis = __o.axis
                    new_ax_id = __o.ax_id
                else:
                    new_axis = []
                    new_ax_id = []

            new_parameter.axis = new_axis
            new_parameter.ax_id = new_ax_id

            return new_parameter

        elif isinstance(__o, numbers.Number):
            new_value = self.value * __o
            new_name = f"{self.name} + {__o}"
            new_parameter = Parameter(
                name=new_name, value=new_value, unit=self.unit)
            if self.axis is not []:
                new_axis = copy.deepcopy(self.axis)
                for i,axis in enumerate(new_axis):
                    new_axis[i]["axis"] = np.array(axis["axis"])* __o
                # new_axis = [np.array([item * __o for item in axis]) for axis in self.axis ]
                new_ax_id = self.ax_id
                new_parameter.axis = new_axis
                new_parameter.ax_id = new_ax_id
            return new_parameter
        
        elif isinstance(__o, np.ndarray):
            if self.axis.shape != __o.shape:
                raise RuntimeError(
                    "Both parameters axis and the array must have the same shape")

    def __rmul__(self, __o:object):
        return self.__mul__(__o)

        





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
        Pulse.__init__(self, t=t, tp=tp,freq=freq, **kwargs)
        self.scale = None
        pass


class Delay(Pulse):
    """
    DEPRECATION WARNING: THIS WILL BE REMOVED SOON.
    
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

