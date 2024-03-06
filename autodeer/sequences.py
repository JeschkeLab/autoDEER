from autodeer.classes import Parameter
from autodeer.pulses import Pulse, Detection, Delay, RectPulse
import autodeer.pulses as ad_pulses
import numpy as np
import matplotlib.pyplot as plt
import json
from itertools import product
import copy
from autodeer.utils import build_table, autoEPRDecoder
import json
from autodeer import __version__
import uuid
import base64
import autodeer.pulses as ad_pulses
import numbers

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
        self.axes_uuid = []
        self.reduce_uuid = []



        if isinstance(B, Parameter):
            self.B = B.copy()
        else:
            self.B = Parameter(
                "B", B, "Gauss",
                "The static B0 field for the experiment")
        
        self.LO = Parameter(
            "LO", LO, "GHz",
            "The local oscillator frequency.")
        
        if isinstance(reptime, Parameter):
            self.reptime = reptime.copy()
        else:
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
        if hasattr(self,'evo_params'):
            acqs *= np.prod([np.prod(param.dim) for param in self.evo_params])        
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
        self.evo_params = params
        self.axes_uuid = [param.uuid for param in params]
        self.reduce_uuid = [param.uuid for param in reduce]


        # Build the ProgTable
        progTable = {"EventID": [], "Variable": [], "axis": [],
                     "axID": [], "uuid": [], "reduce": []}
        
        for n, pulse in enumerate(self.pulses):
            table = pulse.build_table()
            for i in range(len(table["uuid"])):
                if table["uuid"][i] in self.axes_uuid:
                    progTable["axID"].append(self.axes_uuid.index(table["uuid"][i]))
                    progTable["uuid"].append(table["uuid"][i]) 
                    progTable["EventID"].append(n)
                    progTable["Variable"].append(table["Variable"][i])
                    progTable["axis"].append(table["axis"][i])
                    if table["uuid"][i] in self.reduce_uuid:
                        progTable["reduce"].append(True)
                    else:
                        progTable["reduce"].append(False)

        for var_name in vars(self):
            var = getattr(self, var_name)
            if type(var) is Parameter:
                if not var.is_static() and not var.virtual:
                    for i in range(len(var.axis)):
                        if var.axis[i]["uuid"] in self.axes_uuid:
                            progTable["axID"].append(self.axes_uuid.index(var.axis[i]["uuid"]))
                            progTable["EventID"].append(None)
                            progTable["Variable"].append(var_name) 
                            progTable["axis" ].append(var.axis[i]["axis"])
                            progTable["uuid"].append(var.axis[i]["uuid"]) 
                            if var.axis[i]["uuid"] in self.reduce_uuid:
                                progTable["reduce"].append(True)
                            else:
                                progTable["reduce"].append(False)
        self.progTable = progTable
        self._estimate_time()
        return self.progTable
    
    @property
    def seqtable_steps(self):
        if len(self.evo_params) > 0:
            return self.pcyc_dets.shape[0] * (len(self.pcyc_vars)+1) * np.prod([np.prod(param.dim) for param in self.evo_params])

    def shift_detfreq_to_zero(self):
        det_pulse = None
        for pulse in self.pulses:
            if isinstance(pulse,Detection):
                det_pulse = pulse
        
        det_freq = det_pulse.freq.value
        self.LO.value -= det_freq
        for pulse in self.pulses:
            if hasattr(pulse,'freq'):
                pulse.freq.value -= det_freq
            if hasattr(pulse,'init_freq'):
                pulse.init_freq.value -= det_freq
            if hasattr(pulse,'final_freq'):
                pulse.final_freq.value -= det_freq
        return self
    
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
                if param.unit is None:
                    unit = ""
                else:
                    unit = param.unit
                if type(param.value) is str:
                    seq_param_string += "{:<10} {:<12} {:<10} {:<30} \n".format(
                        param.name, param.value, unit, param.description)
                else:
                    seq_param_string += "{:<10} {:<12.5g} {:<10} {:<30} \n".format(
                        param.name, param.value, unit, param.description)
        
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
    
    def _to_dict(self):
        to_return = {"version": __version__, "type": "Sequence", "subclass": str(type(self))}

        for key, var in vars(self).items():
            if isinstance(var, Parameter):
                to_return[key] = var._to_dict()
            if key == "pulses":
                new_list = []
                for pulse in var:
                    new_list.append(pulse._to_dict())
                to_return[key] = new_list
            else:
                to_return[key] = var

        return to_return

    def _to_json(self):
        class autoEPREncoder(json.JSONEncoder):
            def default(self, obj):
                if isinstance(obj, np.ndarray):
                    if (len(obj) > 0 ) and isinstance(obj[0], str):
                        return list(obj)
                    data = np.ascontiguousarray(obj.data)
                    data_b64 = base64.b64encode(data)
                    return dict(__ndarray__=str(data_b64),
                                dtype=str(obj.dtype),
                                shape=obj.shape)
                if isinstance(obj, complex):
                    return str(obj)
                if isinstance(obj, numbers.Number):
                    return str(obj)
                if isinstance(obj, uuid.UUID):
                    return_dict = {"__uuid__": str(obj)}
                    return return_dict
                if isinstance(obj, Parameter):
                    return obj._to_dict()
                if isinstance(obj, Pulse):
                    return obj._to_dict()
                if isinstance(obj, Sequence):
                    return obj._to_dict()
                else:
                    return json.JSONEncoder.default(self, obj)
        
        return json.dumps(self._to_dict(), cls=autoEPREncoder, indent=4)

    def save(self, filename):
        """Save the sequence to a JSON file.

        Parameters
        ----------
        filename : str
            Path to the JSON file.

        Returns
        -------
        None

        Raises
        ------
        TypeError
            If the object cannot be serialized to JSON.

        Example
        -------
        >>> obj = Sequence()
        >>> obj.save("my_sequence.json")
        """
        
                
        with open(filename, "w") as f:
           f.write(self._to_json())

    @classmethod
    def _from_dict(cls, dct):
        name = dct["name"]
        B = Parameter._from_dict(dct["B"])
        LO = Parameter._from_dict(dct["LO"])
        reptime = Parameter._from_dict(dct["reptime"])
        averages = Parameter._from_dict(dct["averages"])
        shots = Parameter._from_dict(dct["shots"])
        new_sequence = cls(
            name=name, B=B, LO=LO, reptime=reptime, averages=averages, shots=shots
        )
        for key, var in dct.items(): 
            if isinstance(var, dict) and ("type" in var):
                setattr(new_sequence, key, Parameter._from_dict(var))
            elif key == "pulses":
                for pulse in var:
                    if hasattr(ad_pulses,pulse['type']):
                        new_sequence.pulses.append(
                            getattr(ad_pulses,pulse['type'])._from_dict(pulse))
                    else:
                        new_sequence.pulses.append(
                            Pulse._from_dict(pulse))
            else:
                setattr(new_sequence, key, var)

        return new_sequence
    
    @classmethod
    def _from_json(cls, JSONstring):
        dct = json.loads(JSONstring, object_hook=autoEPRDecoder)
        return cls._from_dict(dct)
    
    @classmethod
    def load(cls, filename):
        """Load an object from a JSON file.

        Parameters
        ----------
        filename : str
            Path to the JSON file.

        Returns
        -------
        obj : Sequence
            The Sequence loaded from the JSON file.

        Raises
        ------
        FileNotFoundError
            If the file does not exist.

        Example
        -------
        >>> obj = Sequence.load("my_sequence.json")
        """
        with open(filename, "r") as f:
           file_buffer = f.read()
        return cls._from_json(file_buffer)

# =============================================================================
#                                Subclasses
# =============================================================================

class DEERSequence(Sequence):
    
    """
    Represents a DEER/PELDOR sequence. 
    """

    def __init__(
        self, *, tau1, tau2, tau3 = None, tau4=None, dt, B, LO, reptime,
        averages, shots, ESEEM_avg=None, **kwargs) -> None:
        """Build a DEER sequence using rectangular pulses

        Parameters
        ----------
        tau1 : int or float
            The first interpulse delay in us
        tau2 : int or float
            The second interpulse delay in us
        dt : int or float
            The time step for DEER measurment in ns
        B : int or float
            The B0 field, in Guass
        LO : int or float
            The LO frequency in GHz
        reptime : _type_
            The shot repetition time in us
        averages : int
            The number of scans.
        shots : int
            The number of shots per point
        tau3 : int or float, optional
            The delay between the first static pump pulse in 5-pulse DEER and 
            the 1st refocusing pulse in us, by default None. If the value is 
            None then a 4-pulse sequence is created, otherwise a 5-pulse. 
        ESEEM_avg: str
            Selecting the ESEEM averaging required. ESEEM averaging works by 
            introducing small stepping in the first tau delay and averaging 
            across them all. Options:
            * `proton` - 8 small steps of 8ns
            * `deuteron` - 8 small steps of 16ns
            * None - default
        
        Optional Parameters
        -------------------
        exc_pulse : Pulse
            An autoEPR Pulse object describing the excitation pi/2 pulse. If
            not specified a RectPulse will be created instead. 
        ref_pulse : Pulse
            An autoEPR Pulse object describing the refocusing pi pulses. If
            not specified a RectPulse will be created instead. 
        pump_pulse : Pulse
            An autoEPR Pulse object describing the pumping pi pulse/s. If
            not specified a RectPulse will be created instead. 
        det_event : Pulse
            An autoEPR Pulse object describing the detection even. If
            not specified a standard detection event will be created instead,
            with an offset frequency of 0MHz. 

        """
        name = "DEER sequence"
        self.tau1us = tau1
        self.tau1 = Parameter(name="tau1", value=tau1 * 1e3, unit="ns", description="The first interpulse delays")
        self.tau2 = Parameter(name="tau2", value=tau2 * 1e3, unit="ns", description="The second interpulse delays")
        if tau3 is not None:
            self.tau3 = Parameter(name="tau3", value=tau3 * 1e3, unit="ns", description="The delay between static pump and the fisrt refocusing pulse.")
        if tau4 is not None:
            self.tau4 = tau4 * 1e3

        self.dt = dt
        self.deadtime = 200
        self.ESEEM = False
        self.add_ESEEM_avg(ESEEM_avg)

        super().__init__(
            name=name, B=B, LO=LO, reptime=reptime, averages=averages,
            shots=shots, **kwargs)

        if "exc_pulse" in kwargs:
            self.exc_pulse = kwargs["exc_pulse"]
        if "ref_pulse" in kwargs:
            self.ref_pulse = kwargs["ref_pulse"]
        if "pump_pulse" in kwargs:
            self.pump_pulse = kwargs["pump_pulse"]
        if "det_event" in kwargs:
            self.det_event = kwargs["det_event"]
        
    def add_ESEEM_avg(self, type=None):
        if type is None:
            self.tau1 = Parameter(name="tau1", value=self.tau1us * 1e3, unit="ns", description="The first interpulse delays",virtual=True)
            self.ESEEM=False
        elif (type.lower() == 'proton') or (type.lower() == 'p') or (type.lower() == 'h'):
            self.tau1 = Parameter(
                name="tau1", value=self.tau1us * 1e3, unit="ns", step=8, dim=8,
                description="The first interpulse delays with proton ESEEM averaging",virtual=True)
            self.ESEEM=True
        elif (type.lower() == 'deuteron') or (type.lower() == 'd'):
            self.tau1 = Parameter(
                name="tau1", value=self.tau1us * 1e3, unit="ns", step=16, dim=8,
                description="The first interpulse delays with deuteron ESEEM averaging",virtual=True)
            self.ESEEM=True


    def three_pulse(self, tp=16):
        """
        Build a four pulse DEER sequence.

        Parameters
        ----------
        tp : float
            Length of default RectPulse in ns, by default 16ns.
        """

        self.name = "3pDEER"

        dim = np.floor((self.tau1 - 2*self.deadtime)/self.dt)

        t = Parameter(name="t", value=self.deadtime, step=self.dt, dim=dim, unit="ns", description="The time axis")

        if hasattr(self,"exc_pulse"): # Exc pulse
            self.addPulse(self.exc_pulse.copy(t=0))
        else:
            self.addPulse(RectPulse(  
                t=0, tp=tp, freq=0, flipangle=np.pi/2
            ))
        

        if hasattr(self,"pump_pulse"): # Pump 1 pulse
            self.addPulse(self.pump_pulse.copy(t=t))
        else:
            self.addPulse(RectPulse(
                t=t, tp=tp, freq=0, flipangle=np.pi
            ))

        if hasattr(self,"ref_pulse"): # Ref 2 pulse
            self.addPulse(self.ref_pulse.copy(t=self.tau1))
        else:
            self.addPulse(RectPulse( 
                t=self.tau1, tp=tp, freq=0, flipangle=np.pi
            ))

        if hasattr(self, "det_event"):
            self.addPulse(self.det_event.copy(t=2*self.tau1))
        else:
            self.addPulse(Detection(t=2*self.tau1, tp=512))


        self.evolution([t])

        # if self.ESEEM_avg.lower() == "proton":
        #     ESEEM_axis = np.arange(0,8) * 8
        #     self.addPulsesProg(
        #         [2, 3],
        #         ["t"],
        #         1,
        #         [self.tau1+ ESEEM_axis, 2*(self.tau1+ ESEEM_axis)])
        # elif self.ESEEM_avg.lower() == "deuteron":
        #     ESEEM_axis = np.arange(0,8) * 16
        #     self.addPulsesProg(
        #         [2, 3],
        #         ["t"],
        #         1,
        #         [self.tau1+ ESEEM_axis, 2*(self.tau1+ ESEEM_axis)])

        pass

    def four_pulse(self, tp=16, relaxation=False):
        """
        Build a four pulse DEER sequence.

        Parameters
        ----------
        tp : float
            Length of default RectPulse in ns, by default 16ns.
        """
        self.name = "4pDEER"
        if relaxation:
            self.tau1 = Parameter(
                name="tau1", value=self.tau1.value, unit="ns",
                description="The first interpulse delays", virtual=True)
            self.tau2 = Parameter(
                name="tau2", value=self.tau2.value, unit="ns", dim=100,
                step=50,
                description="The second interpulse delays", virtual=True)
            self.t = Parameter(name="t", value=0, unit="ns", description="The time axis", virtual=True)
        else:
            dim = np.floor((self.tau2.value)/self.dt)
            self.t = Parameter(name="t", value = - self.deadtime, step=self.dt,
                       dim=dim, unit="ns", description="The time axis", virtual=True)

        
        if hasattr(self,"exc_pulse"): # Exc pulse
            self.addPulse(self.exc_pulse.copy(t=0))
        else:
            self.addPulse(RectPulse(  
                t=0, tp=tp, freq=0, flipangle=np.pi/2
            ))
        
        
        if hasattr(self,"ref_pulse"): # Ref 1 pulse
            self.addPulse(self.ref_pulse.copy(t=self.tau1))
        else:
            self.addPulse(RectPulse(t=self.tau1, tp=tp, freq=0, flipangle=np.pi))
        
        if hasattr(self,"pump_pulse"): # Pump 1 pulse
            self.addPulse(self.pump_pulse.copy(t=2*self.tau1 + self.t))
        else:
            self.addPulse(RectPulse(t=2*self.tau1 + self.t,
                                    tp=tp, freq=0, flipangle=np.pi))


        if hasattr(self,"ref_pulse"): # Ref 2 pulse
            
            self.addPulse(self.ref_pulse.copy(t=2*self.tau1 + self.tau2))
        else:
            self.addPulse(RectPulse(t=2*self.tau1 + self.tau2,
                                    tp=tp, freq=0, flipangle=np.pi))

        if hasattr(self, "det_event"):
            self.addPulse(self.det_event.copy(t=2*(self.tau1+self.tau2)))
        else:
            self.addPulse(Detection(t=2*(self.tau1+self.tau2), tp=512))


        if relaxation:
            self.evolution([self.tau2])
        elif self.ESEEM:
            self.evolution([self.t,self.tau1],[self.tau1])
        else:
            self.evolution([self.t])


    def five_pulse(self, tp=16, relaxation=False, re_step=50, re_dim=100):
        """
        Build a five pulse DEER sequence.

        Parameters
        ----------
        tp : float
            Length of default RectPulse in ns, by default 16ns.
        """
        

        # t = Parameter(name="t", value=-self.deadtime, step=self.dt,
        #                dim=dim, unit="ns", description="The time axis")
        
        self.name = "5pDEER"
        if relaxation:
            self.tau1 = Parameter(
                name="tau1", value=self.tau1.value, unit="ns", dim=re_dim, step=re_step,
                description="The first interpulse delays", virtual=True)
            self.tau2 = Parameter(
                name="tau2", value=self.tau2.value, unit="ns", dim=re_dim,
                step=re_step, link=self.tau1,
                description="The second interpulse delays", virtual=True)
            
            self.t = Parameter(
                name="t", value=self.tau3.value, unit="ns", description="The time axis",
                virtual=True)

        else:
            if self.deadtime > self.tau3.value:
                raise ValueError("Deadtime must be greater than tau3, otherwise pulses crash")
            
            dim = np.floor((self.tau2.value + self.tau1.value - self.tau3.value)/self.dt)
            self.t = Parameter(
                name="t", value= self.tau3.value - self.deadtime, step=self.dt,
                dim=dim, unit="ns", description="The time axis", virtual=True)




        if hasattr(self,"exc_pulse"): # Exc pulse
            self.addPulse(self.exc_pulse.copy(t=0))
        else:
            self.addPulse(RectPulse(  
                t=0, tp=tp, freq=0, flipangle=np.pi/2
            ))

        if hasattr(self,"pump_pulse"): # Pump 1 pulse
            self.addPulse(self.pump_pulse.copy(t=self.tau1 - self.tau3))
        else:
            self.addPulse(RectPulse(t=self.tau1 - self.tau3,
                                    tp=tp, freq=0, flipangle=np.pi))

        if hasattr(self,"ref_pulse"): # Ref 1 pulse
            self.addPulse(self.ref_pulse.copy(t=self.tau1))
        else:
            self.addPulse(RectPulse(t=self.tau1, tp=tp, freq=0,
                                    flipangle=np.pi))
        
        if hasattr(self,"pump_pulse"): # Pump 2 pulse
            self.addPulse(self.pump_pulse.copy(t=self.tau1 + self.t))
        else:
            self.addPulse(RectPulse(t=self.tau1 + self.t, tp=tp, freq=0,
                                    flipangle=np.pi))

        if hasattr(self,"ref_pulse"): # Ref 2 pulse
            self.addPulse(self.ref_pulse.copy(t=2*self.tau1 + self.tau2))
        else:
            self.addPulse(RectPulse(t=2*self.tau1 + self.tau2, tp=tp, freq=0,
                                    flipangle=np.pi))

        if hasattr(self, "det_event"):
            self.addPulse(self.det_event.copy(t=2*(self.tau1+self.tau2)))
        else:
            self.addPulse(Detection(t=2*(self.tau1+self.tau2), tp=512))



        if relaxation:
            self.evolution([self.tau1])

        else:
            self.evolution([self.t])

        # self.addPulsesProg(
        #     [3],
        #     ["t"],
        #     0,
        #     axis)
        
        # if self.ESEEM_avg is not None:
        #     # This is not quite a perfect ESEEM average. I am not sure if this
        #     # is possible with a pulse based scripting. I think only following
        #     # a delay based model can you easily describe two axes acting on 
        #     # the same pulse. 
        #     if self.ESEEM_avg.lower() == "proton":
        #         ESEEM_axis = np.arange(0,8) * 8
        #     elif self.ESEEM_avg.lower() == "deuteron":
        #         ESEEM_axis = np.arange(0,8) * 16
        #     self.addPulsesProg(
        #         [1, 2, 4, 5],
        #         ["t"],
        #         1,
        #         [p1 - ESEEM_axis, r1 + 1*ESEEM_axis, r2 + 2*ESEEM_axis,
        #          d + 2*ESEEM_axis])

        # pass

    def seven_pulse(self, tp=16, relaxation=False):
        """
        Build a seven pulse DEER sequence.

        Parameters
        ----------
        tp : float
            Length of default RectPulse in ns, by default 16ns.
        """
        self.name = "7pDEER"

        if relaxation:
            self.tau1 = Parameter(
                name="tau1", value=self.tau1.value, unit="ns", dim=100, step=50,
                description="The first interpulse delays", virtual=True)
            self.tau2 = Parameter(
                name="tau2", value=self.tau2.value, unit="ns", dim=100,
                step=50, link=self.tau1,
                description="The second interpulse delays", virtual=True)
            
            self.t = Parameter(
                name="t", value=0, unit="ns", description="The time axis",
                virtual=True)

        else:
            dim = np.floor((self.tau2.value + self.tau1.value - 2*self.deadtime)/self.dt)
            self.t = Parameter(
                name="t", value = + self.deadtime, step=self.dt,
                dim=dim, unit="ns", description="The time axis", virtual=True)

        # axis = np.arange(0, self.tau2 + self.tau1 - 2*self.deadtime, self.dt)
        
        # dim = np.floor((self.tau2 + self.tau1 - 2*self.deadtime)/self.dt)

        # t = Parameter(name="t", value=0, step=self.dt,
        #                dim=dim, unit="ns", description="The time axis")


        if hasattr(self,"exc_pulse"): # Exc pulse
            self.addPulse(self.exc_pulse.copy(t=0))
        else:
            self.addPulse(RectPulse(  
                t=0, tp=tp, freq=0, flipangle=np.pi/2
            ))

        if hasattr(self,"ref_pulse"): # Ref 1 pulse
            self.addPulse(self.ref_pulse.copy(t=self.tau1))
        else:
            self.addPulse(RectPulse(t=self.tau1, tp=tp, freq=0, flipangle=np.pi))

        if hasattr(self,"pump_pulse"): # Pump 1 pulse
            self.addPulse(self.pump_pulse.copy(t=self.tau1 + self.t))
        else:
            self.addPulse(RectPulse(t=self.tau1 + self.t, tp=tp, freq=0, flipangle=np.pi))

        if hasattr(self,"ref_pulse"): # Ref 2 pulse
            self.addPulse(self.ref_pulse.copy(t= 2*self.tau1 + self.tau2))
        else:
            self.addPulse(RectPulse( 
                t= 2*self.tau1 + self.tau2, tp=tp, freq=0, flipangle=np.pi
            ))

        if hasattr(self,"pump_pulse"): # Pump 2 pulse
            self.addPulse(self.pump_pulse.copy(t=2*self.tau1 + self.tau2 + self.tau4))
        else:
            self.addPulse(RectPulse(t=2*self.tau1 + self.tau2 + self.tau4, tp=tp, freq=0, flipangle=np.pi))
        
        r3_pos = 2*self.tau1 + 2*self.tau2 + self.tau3
        if hasattr(self,"pump_pulse"): # Pump 3 pulse
            self.addPulse(self.pump_pulse.copy(t=r3_pos - self.t))
        else:
            self.addPulse(RectPulse(t=r3_pos -self.t, tp=tp, freq=0, flipangle=np.pi))

        
        if hasattr(self,"ref_pulse"): # Ref 3 pulse
            self.addPulse(self.ref_pulse.copy(t=r3_pos))
        else:
            self.addPulse(RectPulse( 
                t=r3_pos, tp=tp, freq=0, flipangle=np.pi))

        d_pos = 2*(self.tau1+self.tau2+self.tau3)
        if hasattr(self, "det_event"):
            self.addPulse(self.det_event.copy(t=d_pos))
        else:
            self.addPulse(Detection(t=d_pos, tp=512))


        if relaxation:
            self.evolution([self.tau1])

        else:
            self.evolution([self.t])

        # self.addPulsesProg(
        #     pulses=[2],
        #     variables=["t"],
        #     axis_id=0,
        #     axis=axis+self.pulses[2].t.value)
        
        # self.addPulsesProg(
        #     pulses=[5],
        #     variables=["t"],
        #     axis_id=0,
        #     axis=axis+self.pulses[2].t.value)
        pass
    
    def nDEER_CP(self, n:int, tp=16, relaxation=False, pcyc='Normal'):
        """Generate an nDEER sequence.

        The sum of tau1 and tau2 is used as total trace length. 


        Parameters
        ----------
        n : int
            The number of refocusing pulses
        tp : int, optional
            _description_, by default 16
        relaxation : bool, optional
            _description_, by default False
        pcyc: str, optional
            Normal: Phases cycles pump and observer pulses, no DC cycle
            NormalDC: Phases cycles pump and observer pulses, DC cycle
            Obs: Phases cycles observer pulses, no DC cycle
            ObsDC: Phases cycles and observer pulses, DC cycle
        """

        self.name = "nDEER-CP"
        self.n = Parameter(
            name="n", value=n,
            unit=None, description="The number of refoucsing pulses",
            virtual=True)

        self.step = Parameter(
            name="step", value=2*(self.tau1.value+self.tau2.value)/(2*n),
            unit="ns", description="The gap of the first refoucsing pulse",
            virtual=True)

        if relaxation:

            self.step = Parameter(
                name="step", value=300,
                dim=100, step=200,
                unit="ns", description="The gap of the first refoucsing pulse",
                virtual=True)            
            self.t = Parameter(
                name="t", value=0, unit="ns", description="The time axis",
                virtual=True)

        else:
            
            dim = np.floor((n*self.step.value)/self.dt)

            self.t = Parameter(
                name="t", value= - self.deadtime, step=self.dt,
                dim=dim, unit="ns", description="The time axis", virtual=True)
                

        if hasattr(self,"exc_pulse"): # Exc pulse
            self.addPulse(self.exc_pulse.copy(t=0))
        else:
            self.addPulse(RectPulse(  
                t=0, tp=tp, freq=0, flipangle=np.pi/2
            ))
        
        if hasattr(self,"ref_pulse"): # Ref 1 pulse
            self.addPulse(self.ref_pulse.copy(t=self.step))
        else:
            self.addPulse(RectPulse( 
                t=self.step, tp=tp, freq=0, flipangle=np.pi
            ))

        if hasattr(self,"pump_pulse"): # Pump 1 pulse
            self.addPulse(self.pump_pulse.copy(t=n*self.step + self.t))
        else:
            self.addPulse(RectPulse(
                t=n*self.step + self.t, tp=tp, freq=0, flipangle=np.pi
            ))
        
        for i in np.arange(1,n):

            if hasattr(self,"ref_pulse"): # Ref i pulse
                self.addPulse(self.ref_pulse.copy(t=(2*i+1)*self.step))
            else:
                self.addPulse(RectPulse( 
                    t=(2*i+1)*self.step, tp=tp, freq=0, flipangle=np.pi
                ))

        if hasattr(self, "det_event"):
            self.addPulse(self.det_event.copy(t=2*n*(self.step)))
        else:
            self.addPulse(Detection(t=2*n*(self.step), tp=512))
        
        if relaxation:
            self.evolution([self.step])
        else:
            self.evolution([self.t])
      
        if pcyc == 'NormalDC' or pcyc=='ObsDC':
            self.pulses[0]._addPhaseCycle(
                [0, np.pi], [1, -1])
        
        if pcyc == 'NormalDC' or pcyc == 'Normal':
            self.pulses[2]._addPhaseCycle(
                [0, np.pi/2, np.pi, -np.pi/2], [1, 1, 1, 1])
        
        if pcyc is not None:
            if n-1 == 1:
                self.pulses[1]._addPhaseCycle(
                    [0, np.pi], [1,1])
            elif n > 2:
                self.pulses[1]._addPhaseCycle(
                    [0, np.pi/2, np.pi, -np.pi/2], [1, -1, 1, -1])
                self.pulses[n]._addPhaseCycle(
                    [0, np.pi], [1,1])
                for i in range(3,n):
                    self.pulses[i]._addPhaseCycle(
                    [0, np.pi/2, np.pi, -np.pi/2], [1, -1, 1, -1])

        pass
    def select_pcyc(self, option: str):
        """Choose which DEER phase you would like.
        
        .. |xp| replace:: x\ :sub:`p`

        .. table::
            :width: 150
            :widths: 10 10 10 5 30 30 5
 
            +---------------------------+-------------+----------------+--------+---------------------------+-----------------------------+------------+
            | Phase cycle               | Short Code  | Sequence       | Steps  | Pulse Phase Cycle         | Remaining Echoes            | Ref.       |
            +===========================+=============+================+========+===========================+=============================+============+
            | (x)x|xp|x                 | DC          | ALL            | 2      | [+(+x)-(-x)]              | PE12rp, SE(PE12)p3, PE12rpr3|            |
            +---------------------------+-------------+----------------+--------+---------------------------+-----------------------------+------------+
            | (x)[|xp|]x                | 8step_3p    | 3 pulse        | 8      | [+(+x)-(-x)]              |                             |            |
            +---------------------------+-------------+----------------+--------+---------------------------+-----------------------------+------------+
            | x[x][|xp|]x               | 16step_4p   | 4 pulse        | 16     | [+(+x)-(+y)+(-x)-(-y)]    |                             | [1]        |
            +---------------------------+-------------+----------------+--------+                           +-----------------------------+------------+
            | x|xp|[x][|xp|]x           | 16step_5p   | 5 pulse        | 16     | [+(+x)+(+y)+(-x)+(-y)]    | PEp02r3,b PE1p0r2r3b        | [1]        |
            +---------------------------+-------------+----------------+--------+---------------------------+-----------------------------+------------+
            | x[x]|xp|(x)(|xp|)(|xp|)x  | 32step_7p   | 7 pulse        | 32     |                           |                             | [1]        |
            +---------------------------+-------------+----------------+--------+---------------------------+-----------------------------+------------+

        Parameters
        ----------
        option : str
            The short code of the phase cycle. See table above.
        """

        if option == "DC":
            self.pulses[0]._addPhaseCycle([0,np.pi],[1,-1])


        elif option == "8step_3p":
            self.pulses[0]._addPhaseCycle([0],[1])
            self.pulses[1]._addPhaseCycle(
                [0, np.pi/2, np.pi, -np.pi/2], [1, 1, 1, 1])
        
        elif option == "16step_4p":
            self.pulses[0]._addPhaseCycle([0],[1])
            self.pulses[1]._addPhaseCycle(
                [0, np.pi/2, np.pi, -np.pi/2], [1, -1, 1, -1])
            self.pulses[2]._addPhaseCycle(
                [0, np.pi/2, np.pi, -np.pi/2], [1, 1, 1, 1])
            self.pulses[3]._addPhaseCycle([0],[1])

            
        elif option == "16step_5p":
            self.pulses[0]._addPhaseCycle([0],[1])
            self.pulses[1]._addPhaseCycle([0],[1])
            self.pulses[2]._addPhaseCycle(
                [0, np.pi/2, np.pi, -np.pi/2], [1, -1, 1, -1])
            self.pulses[3]._addPhaseCycle(
                [0, np.pi/2, np.pi, -np.pi/2], [1, 1, 1, 1])
            self.pulses[4]._addPhaseCycle([0],[1])

        elif option =="32step_7p":
            self.pulses[0]._addPhaseCycle([0],[1])                  # Exc
            self.pulses[1]._addPhaseCycle(                          # Ref 1
                [0, np.pi/2, np.pi, -np.pi/2], [1, -1, 1, -1])
            self.pulses[2]._addPhaseCycle([0],[1])                  # Pump 1
            self.pulses[3]._addPhaseCycle([0, np.pi], [1, 1])       # Ref 2
            self.pulses[4]._addPhaseCycle([0, np.pi], [1, 1])       # Pump 2
            self.pulses[5]._addPhaseCycle([0, np.pi], [1, 1])       # Pump 3
            self.pulses[6]._addPhaseCycle([0],[1])                  # Ref 3

        
        self._buildPhaseCycle()


    def simulate(self):
        t = np.arange(0, self.tau1+self.tau2, self.dt)
        r = np.arange(1.5, 8)
        pass

# =============================================================================

class HahnEchoSequence(Sequence):
    """
    Represents a Hahn-Echo sequence. 
    """
    def __init__(self, *, B, LO, reptime, averages, shots, **kwargs) -> None:
        """Build a Hahn-Echo sequence using either rectangular pulses or
        specified pulses. By default no progression is added to this sequence.

        Parameters
        ----------
        B : int or float
            The B0 field, in Guass
        LO : int or float
            The LO frequency in GHz
        reptime : _type_
            The shot repetition time in us
        averages : int
            The number of scans.
        shots : int
            The number of shots per point
        
        Optional Parameters
        -------------------
        pi2_pulse : Pulse
            An autoEPR Pulse object describing the excitation pi/2 pulse. If
            not specified a RectPulse will be created instead. 
        pi_pulse : Pulse
            An autoEPR Pulse object describing the refocusing pi pulses. If
            not specified a RectPulse will be created instead. 
        det_event : Pulse
            An autoEPR Pulse object describing the detection even. If
            not specified a standard detection event will be created instead,
            with an offset frequency of 0MHz. 

        """
        
        name = "Hahn Echo"
        
        if "pi_pulse" in kwargs:
            self.pi_pulse = kwargs["pi_pulse"]
        if "pi2_pulse" in kwargs:
            self.pi2_pulse = kwargs["pi2_pulse"]
        if "det_event" in kwargs:
            self.det_event = kwargs["det_event"]


        super().__init__(
            name=name, B=B, LO=LO, reptime=reptime, averages=averages,
            shots=shots, **kwargs)

        if "tau" in kwargs:
            tau = kwargs["tau"]
        else:
            tau = 500
        if "tp" in kwargs:
            tp = kwargs["tp"]
        else:
            tp = 12

        if hasattr(self, "pi2_pulse"):
            self.addPulse(self.pi2_pulse.copy(
                t=0, pcyc={"phases":[0, np.pi], "dets": [1, -1]}))
        else:
            self.addPulse(RectPulse(  # Exc pulse
                t=0, tp=tp, freq=0, flipangle=np.pi/2, 
                pcyc={"phases":[0, np.pi], "dets":[1, -1]}
            ))

        if hasattr(self, "pi_pulse"):
            pi_pulse = self.addPulse(self.pi_pulse.copy(
                t=tau, pcyc={"phases":[0], "dets": [1]}))
        else:
            pi_pulse = self.addPulse(RectPulse( # Pump 1 pulse
                t=tau, tp=tp, freq=0, flipangle=np.pi
            ))

        if hasattr(self, "det_event"):
            self.addPulse(self.det_event.copy(t=2*tau))
        else:
            self.addPulse(Detection(t=2*tau, tp=self.det_window.value))

# =============================================================================

class T2RelaxationSequence(HahnEchoSequence):
    """
    Represents a T2 relaxation sequence. A Hahn Echo where the interpulse delay increases
    """

    def __init__(self, *, B, LO, reptime, averages, shots, step=40, dim=200, **kwargs) -> None:

        self.tau = Parameter(name="tau", value=500,step=step,dim=dim, unit="ns", description="The interpulse delay",virtual=True)
        super().__init__(B=B, LO=LO, reptime=reptime, averages=averages, shots=shots,tau=self.tau, **kwargs)

        self.name = "T2 Relaxation"
        self.evolution([self.tau])

# =============================================================================

class FieldSweepSequence(HahnEchoSequence):
    """
    Represents a Field Sweep (EDFS) sequence. 
    """
    def __init__(self, *, B, LO, Bwidth, reptime, averages, shots, **kwargs) -> None:
        """Build a Field Sweep (EDFS) sequence using either rectangular pulses or
        specified pulses.

        Parameters
        ----------
        B : int or float
            The B0 field, in Guass
        Bwidth: int or float
            The width of the field sweep, in Gauss
        LO : int or float
            The LO frequency in GHz
        reptime : _type_
            The shot repetition time in us
        averages : int
            The number of scans.
        shots : int
            The number of shots per point
        
        Optional Parameters
        -------------------
        pi2_pulse : Pulse
            An autoEPR Pulse object describing the excitation pi/2 pulse. If
            not specified a RectPulse will be created instead. 
        pi_pulse : Pulse
            An autoEPR Pulse object describing the refocusing pi pulses. If
            not specified a RectPulse will be created instead. 
        """
        super().__init__(
            B=B, LO=LO, reptime=reptime, averages=averages,
            shots=shots, **kwargs)
        self.name = "Field Sweep"


        self.B = Parameter(
            "B", value=B, start = -Bwidth/2, step=1, dim=Bwidth, unit="Gauss", description="Field sweep width"
        )
        
        self.evolution([self.B])
        
# =============================================================================

class ReptimeScan(HahnEchoSequence):
    """
    Represents a reptime scan of a Hahn Echo Sequence. 
    """
    def __init__(self, *, B, LO, reptime, reptime_max, averages, shots, **kwargs) -> None:
        """A Hahn echo sequence is perfomed with the shot repetition time increasing.1

        Parameters
        ----------
        B : int or float
            The B0 field, in Guass
        LO : int or float
            The LO frequency in GHz
        reptime: float
            The default reptime, this is used for tuning pulses etc...
        reptime_max : np.ndarray
            The maximum shot repetition time in us    
        averages : int
            The number of scans.
        shots : int
            The number of shots per point
        
        Optional Parameters
        -------------------
        pi2_pulse : Pulse
            An autoEPR Pulse object describing the excitation pi/2 pulse. If
            not specified a RectPulse will be created instead. 
        pi_pulse : Pulse
            An autoEPR Pulse object describing the refocusing pi pulses. If
            not specified a RectPulse will be created instead. 
        """
        min_reptime = 20
        dim = 100
        step  = (reptime_max-min_reptime)/dim
        step = np.around(step,decimals=-1)
        step = np.around(step,decimals=-1)
        reptime = Parameter(
            "reptime", reptime, start = min_reptime-reptime, step=step, dim=100, unit="us",
            description = "The shot repetition time")
        
        super().__init__(
            B=B, LO=LO, reptime=reptime, averages=averages,
            shots=shots, **kwargs)
        self.name = "reptime Scan"

        self.evolution([self.reptime])

# =============================================================================

class CarrPurcellSequence(Sequence):
    """
    Represents a Carr-Purcell sequence. 
    """
    def __init__(self, *, B, LO, reptime, averages, shots,
            tau, n, dim=100,**kwargs) -> None:
        """Build a Carr-Purcell dynamical decoupling sequence using either 
        rectangular pulses or specified pulses.

        Parameters
        ----------
        B : int or float
            The B0 field, in Guass
        LO : int or float
            The LO frequency in GHz
        reptime : _type_
            The shot repetition time in us
        averages : int
            The number of scans.
        shots : int
            The number of shots per point
        tau : int
            The maximum total sequence length in us
        n : int
            The number refocusing pulses
        dim : int
            The number of points in the X axis

        Optional Parameters
        -------------------
        pi2_pulse : Pulse
            An autoEPR Pulse object describing the excitation pi/2 pulse. If
            not specified a RectPulse will be created instead. 
        pi_pulse : Pulse
            An autoEPR Pulse object describing the refocusing pi pulses. If
            not specified a RectPulse will be created instead. 
        """

        name = "Carr-Purcell"
        super().__init__(
            name=name, B=B, LO=LO, reptime=reptime, averages=averages,
            shots=shots, **kwargs)
        self.tau = Parameter(name="tau", value=tau*1e3, unit="ns",
            description="Total sequence length", virtual=True)
        self.n = Parameter(name="n", value=n,
            description="The number of pi pulses", unit="None", virtual=True)
        self.dim = Parameter(name="dim", value=dim, unit="None",
            description="The number of points in the X axis", virtual=True)

        if "pi_pulse" in kwargs:
            self.pi_pulse = kwargs["pi_pulse"]
        if "pi2_pulse" in kwargs:
            self.pi2_pulse = kwargs["pi2_pulse"]
        if "det_event" in kwargs:
            self.det_event = kwargs["det_event"]

        self._build_sequence()

    def _build_sequence(self):

        n = self.n.value
        deadtime = 300
        # dt = 20
        # dim = np.floor((self.tau.value/(2*self.n.value) -deadtime)/dt)
        dim = self.dim.value
        tau_interval = self.tau.value/(2*n)
        dt = (tau_interval-deadtime)/dim
        step = Parameter("step",deadtime,unit="ns",step=dt, dim=dim)
        # # multipliers = [1]
        # # multipliers += [1+2*i for i in range(1,n)]
        # # multipliers += [2*n]

        # axis = np.arange(deadtime,tau/(2*n),10)

        if hasattr(self, "pi2_pulse"):
            self.addPulse(self.pi2_pulse.copy(
                t=0, pcyc={"phases":[0, np.pi], "dets": [1, -1]}))
        else:
            self.addPulse(RectPulse(  # pi/2
                t=0, tp=16, freq=0, flipangle=np.pi/2,
                pcyc={"phases":[0, np.pi], "dets": [1, -1]}
            ))

        for i in range(n):
            if i==(n-1):
                phases = [0]
                dets = [1]
            elif i == (n-2):
                phases = [0, np.pi]
                dets = [1, 1]
            else:
                phases = [0, np.pi/2, np.pi, -np.pi/2]
                dets = [1,-1,1,-1]
            if hasattr(self, "pi_pulse"):
                self.addPulse(self.pi_pulse.copy(
                    t=step*(2*i + 1), pcyc={"phases":phases, "dets": dets}))
            else:
                self.addPulse(RectPulse(  # pi
                    t=step*(2*i + 1), tp=32, freq=0, flipangle=np.pi,
                    pcyc={"phases":phases, "dets": dets}
                ))
        if hasattr(self, "det_event"):
            self.addPulse(self.det_event.copy(t=step*(2*n)))
        else:
            self.addPulse(Detection(t=step*(2*n), tp=512))
        
        self.evolution([step])

# =============================================================================

class RefocusedEcho2DSequence(Sequence):

    """
    Represents a 2D Refocused-echo Sequence. 
    """
    def __init__(self, *, B, LO, reptime, averages, shots,
            tau, dim=100, step=50, **kwargs) -> None:
        """Build a 2D Refocused-echo sequence using either 
        rectangular pulses or specified pulses.

        Parameters
        ----------
        B : int or float
            The B0 field, in Guass
        LO : int or float
            The LO frequency in GHz
        reptime : _type_
            The shot repetition time in us
        averages : int
            The number of scans.
        shots : int
            The number of shots per point
        tau : int
            The starting value in ns
        dim: int
            The number of points in both the X and Y axis
        step: float
            The step in ns for both the X and Y axis

        Optional Parameters
        -------------------
        pi2_pulse : Pulse
            An autoEPR Pulse object describing the excitation pi/2 pulse. If
            not specified a RectPulse will be created instead. 
        pi_pulse : Pulse
            An autoEPR Pulse object describing the refocusing pi pulses. If
            not specified a RectPulse will be created instead. 
        """

        name = "Carr-Purcell"
        super().__init__(
            name=name, B=B, LO=LO, reptime=reptime, averages=averages,
            shots=shots, **kwargs)
        self.tau1 = Parameter(name="tau1", value=tau, dim=dim, step=step, unit="ns",
            description="1st interpulse delay")
        self.tau2 = Parameter(name="tau2", value=tau, dim=dim, step=step, unit="ns",
            description="2nd interpulse delay")

        if "pi_pulse" in kwargs:
            self.pi_pulse = kwargs["pi_pulse"]
        if "pi2_pulse" in kwargs:
            self.pi2_pulse = kwargs["pi2_pulse"]
        if "det_event" in kwargs:
            self.det_event = kwargs["det_event"]

        self._build_sequence()

    def _build_sequence(self):
    
        if hasattr(self, "pi2_pulse"):
            self.addPulse(self.pi2_pulse.copy(
                t=0, pcyc={"phases":[0, np.pi], "dets": [1, -1]}))
        else:
            self.addPulse(RectPulse(  # pi/2
                t=0, tp=16, freq=0, flipangle=np.pi/2,
                pcyc={"phases":[0, np.pi], "dets": [1, -1]}
            ))

        if hasattr(self, "pi_pulse"):
            self.addPulse(self.pi_pulse.copy(
                t=self.tau1, pcyc={"phases":[0, np.pi/2, np.pi, -np.pi/2], "dets": [1,-1,1,-1]}))
        else:
            self.addPulse(RectPulse(
                t=self.tau1, tp=32, freq=0, flipangle=np.pi,
                pcyc={"phases":[0, np.pi/2, np.pi, -np.pi/2], "dets": [1,-1,1,-1]}
            ))

        if hasattr(self, "pi_pulse"):
            self.addPulse(self.pi_pulse.copy(
                t=2*self.tau1 + self.tau2, pcyc={"phases":[0, np.pi], "dets": [1,1]}))
        else:
            self.addPulse(RectPulse(
                t=2*self.tau1 + self.tau2, tp=32, freq=0, flipangle=np.pi,
                pcyc={"phases":[0, np.pi], "dets": [1,1]}
            ))
        
        if hasattr(self, "det_event"):
            self.addPulse(self.det_event.copy(t=2*(self.tau1 + self.tau2)))
        else:
            self.addPulse(Detection(t=2*(self.tau1 + self.tau2), tp=512))

        self.evolution([self.tau1, self.tau2])


# =============================================================================

class ResonatorProfileSequence(Sequence):
    """
    Builds nutation based Resonator Profile sequence. 
    """

    def __init__(self, *, B, LO, reptime, averages, shots, fwidth=0.3, **kwargs) -> None:
        """Build a resonator profile nutation sequence using either 
        rectangular pulses or specified pulses.

        Parameters
        ----------
        B : int or float
            The B0 field, in Guass
        Bwidth: int or float
            The width of the field sweep, in Gauss
        LO : int or float
            The LO frequency in GHz
        reptime : _type_
            The shot repetition time in us
        averages : int
            The number of scans.
        shots : int
            The number of shots per point
        fwidth: float
            The frequency width of the resonator profile in GHz, 
            by default 0.3GHz
        tau1: float
            The delay between the nutating pulse and the Hahn Echo, 
            by default 2000 ns
        tau2: float
            The interpulse delay in the Hahn Echo, 
            by default 500 ns

        Optional Parameters
        -------------------
        pi2_pulse : Pulse
            An autoEPR Pulse object describing the excitation pi/2 pulse. If
            not specified a RectPulse will be created instead. 
        pi_pulse : Pulse
            An autoEPR Pulse object describing the refocusing pi pulses. If
            not specified a RectPulse will be created instead. 
        """

        name = "Resonator-Profile"
        super().__init__(
            name=name, B=B, LO=LO, reptime=reptime, averages=averages,
            shots=shots, **kwargs)
        self.gyro = LO/B

        if "pi_pulse" in kwargs:
            self.pi_pulse = kwargs["pi_pulse"]
        if "pi2_pulse" in kwargs:
            self.pi2_pulse = kwargs["pi2_pulse"]

        self._build_sequence()

    def _build_sequence(self):

        tau1=2000
        tau2=500

        tp = Parameter("tp", 0, step=2, dim=40, unit="ns", description="Test Pulse length")
        fwidth= 0.3
        fstep = 0.02
        dim = np.floor(fwidth*2/0.02)
        center_LO = self.LO.value
        self.LO = Parameter("LO", center_LO, start=-fwidth, step=fstep, dim=dim, unit="GHz", description="LO frequency")
        self.B = Parameter(
            "B",((self.LO.value)/self.gyro), start=-fwidth/self.gyro, step=fstep/self.gyro, dim=dim,
            unit="Guass",link=self.LO,description="B0 Field" )
        
        self.addPulse(RectPulse(  # Hard pulse
            t=0, tp=tp, freq=0, flipangle="Hard"
        ))

        if hasattr(self, "pi2_pulse"):
            self.addPulse(self.pi2_pulse.copy(
                t=tau1, pcyc={"phases":[0, np.pi], "dets": [1, -1]}))
        else:
            self.addPulse(RectPulse(  # pi/2
            t=tau1, tp=16, freq=0, flipangle=np.pi/2, 
            pcyc={"phases":[0, np.pi], "dets": [1, -1]}
            ))
        
        if hasattr(self, "pi_pulse"):
            self.addPulse(self.pi_pulse.copy(
                t=tau1+tau2))
        else:
            self.addPulse(RectPulse(  # pi/2
            t=tau1+tau2, tp=32, freq=0, flipangle=np.pi
            ))

        self.addPulse(Detection(t=tau1+2*tau2, tp=512))


        self.pulses[0].scale.value = 1
        # nut_axis = np.arange(0,66,2,)
        # self.addPulsesProg(
        #     [0],
        #     ["tp"],
        #     0,
        #     nut_axis)

        # # Add frequency sweep
        # width= 0.3
        # axis = np.arange(self.LO.value-width,self.LO.value+width+0.02,0.02)
        # self.addPulsesProg(
        #     [None, None],
        #     ["LO", "B"],
        #     1,
        #     axis,
        #     multipliers=[1,1/self.gyro])
        
        self.evolution([tp, self.LO])

# =============================================================================

class TWTProfileSequence(Sequence):
    """
    Builds TWT based Resonator Profile sequence. 
    """
    
    def __init__(self,*,B,LO,reptime,averages=1,shots=100,**kwargs) -> None:

        name = "TWT-Profile"
        super().__init__(
            name=name, B=B, LO=LO, reptime=reptime, averages=averages,
            shots=shots, **kwargs)

        self._build_sequence()

    def _build_sequence(self,):

        tau1=2000
        tau2=500

        self.addPulse(RectPulse(  # Hard pulse
            t=0, tp=4, freq=0, flipangle="Hard"
        ))

        self.addPulse(RectPulse(  # pi/2
            t=tau1, tp=16, freq=0, flipangle=np.pi/2
        ))
        self.addPulse(RectPulse(  # pi
            t=tau1+tau2, tp=32, freq=0, flipangle=np.pi
        ))

        self.addPulse(Detection(t=tau1+2*tau2, tp=512))


        self.pulses[0].scale.value = 1
        nut_axis = np.arange(0,66,2)
        self.addPulsesProg(
            [0],
            ["tp"],
            0,
            nut_axis)

        # Add amplitude sweep
        width= 0.3
        axis = np.arange(0,1.01,0.01)
        self.addPulsesProg(
                [0],
                ["scale"],
                1,
                axis)

# =============================================================================
