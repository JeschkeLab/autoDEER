from autodeer.classes import Parameter
from autodeer.utils import build_table, sop, autoEPRDecoder
from memoization import cached
import numpy as np
from scipy.integrate import cumulative_trapezoid
import numbers
import uuid
import json
import base64
import matplotlib.pyplot as plt
import scipy.fft as fft
from autodeer import __version__
import copy


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
        elif scale is None:
            self.scale = Parameter("scale", None, None, "Amplitude of pulse")


        self.name = name
        self.Progression = False

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

    @property
    def bandwidth(self):
        BW = Parameter("Bandwidth", 0, "GHz", "Bandwidth of pulse")
        if hasattr(self, "freq"):
            BW.value =0
        elif hasattr(self, "init_freq") & hasattr(self, "final_freq"):
            BW.value = np.abs(self.final_freq.value - self.init_freq.value)
        else:
            raise ValueError()
        return BW
    
    def _addPhaseCycle(self, phases, detections=None):
            """
            Adds a phase cycle to the pulse sequence.

            Args:
                phases (list): List of phases to add to the phase cycle.
                detections (list, optional): List of detection signs. Defaults to None. If None then all cycles are summed.

            Returns:
                None
            """
            if detections is None:
                detections = np.ones(len(phases))
            self.pcyc = {"Phases": list(phases), "DetSigns": list(detections)}
            pass

    def _buildFMAM(self, func, ax=None):
        """
        Builds the amplitude modulation (AM) and frequency modulation (FM) of a given function.

        Args:
            func: A function that takes in an array of values and returns two arrays, representing the AM and FM of the function.

        Returns:
            Two arrays representing the AM and FM of the function.
        """
        self.ax = np.linspace(-self.tp.value/2, self.tp.value/2, 1001)
        dt = self.ax[1]-self.ax[0]
        self.AM, self.FM = func(self.ax)
        FM_arg = 2*np.pi*cumulative_trapezoid(self.FM, initial=0) * dt
        self.complex = self.AM * (np.cos(FM_arg) +1j* np.sin(FM_arg))
        self.FM
        self.AM
        return self.AM, self.FM
    
    
    def build_shape(self,ax=None):
        if ax is None:
            ax = self.ax
        dt = ax[1]-ax[0]
        pulse_points = int(np.around(self.tp.value/dt))
        total_points = ax.shape[0]
        remaining_points = total_points-pulse_points
        pulse_axis = np.zeros(pulse_points)
        AM, FM = self.func(pulse_axis)
        AM = np.pad(AM,(int(np.floor(remaining_points/2)),int(np.ceil(remaining_points/2))),'constant',constant_values=0)
        FM = np.pad(FM,(int(np.floor(remaining_points/2)),int(np.ceil(remaining_points/2))),'constant',constant_values=0)
        FM_arg = 2*np.pi*cumulative_trapezoid(FM, initial=0) * dt
        shape = AM * (np.cos(FM_arg) +1j* np.sin(FM_arg))
 
        return shape

    def build_table(self):
            """
            Builds a table of variables, axes, and UUIDs for all non-static Parameters in the object.

            Returns:
                dict: A dictionary containing the following keys: "Variable", "axis", and "uuid".
                The values for each key are lists of the corresponding values for each non-static Parameter.
            """
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
            """
            Check if all parameters in the pulse object are static.

            Returns:
                bool: True if all parameters are static, False otherwise.
            """
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
    
    @property
    def amp_factor(self):
        """ The B1 amplitude factor (nutation frequency) for the pulse in GHz"""
        amp_factor_value=  self.flipangle.value / (2 * np.pi * np.trapz(self.AM,self.ax))
        return Parameter("amp_factor", amp_factor_value, "GHz", "Amplitude factor for the pulse")

    # @cached(thread_safe=False)
    def exciteprofile(self, freqs=None, resonator = None):
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
        resonator: ad.ResonatorProfile, optional

        
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

        self._buildFMAM(self.func)

        if freqs is None:
            fxs, fs = self._calc_fft()
            fs = np.abs(fs)
            points = np.argwhere(fs>(fs.max()*0.5))[:,0]
            BW = fxs[points[-1]] - fxs[points[0]]
            c_freq = np.mean([fxs[points[-1]], fxs[points[0]]])
            offsets = np.linspace(-2*BW, 2*BW, 100) + c_freq
        else:
            offsets = freqs

        offsets 
        t = self.ax 
        nOffsets = offsets.shape[0]

        ISignal = np.real(self.complex) * self.amp_factor.value
        QSignal = np.imag(self.complex) * self.amp_factor.value

        if resonator is not None:
            FM = self.FM
            amp_factor = np.interp(FM, resonator.freqs-resonator.LO_c, resonator.profile)
            amp_factor = np.min([amp_factor,np.ones_like(amp_factor)*self.amp_factor.value],axis=0)
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
            Mag[2, iOffset] = -2 * (Sz * Density.T).sum().real
        
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
        param_table =  f"{'Name':^12} \t {'Value':^12} \t {'Unit':^7} \t Description\n"
        param_table += f"{'----':^12} \t {'-----':^12} \t {'----':^7} \t -----------\n"
        for attr_name in dir(self):
            attr = getattr(self, attr_name)
            if type(attr) is Parameter:
                if attr.name == "flipangle":
                    if attr.value == np.pi:
                        value = f"{'π':>12}"
                    elif attr.value == np.pi/2:
                        value = f"{'π/2':>12}"
                    else:
                        value = f"{attr.value:>5.5g}"
                elif attr.value is None:
                    value = f"{'N/A':>12}"
                else:
                    value = f"{attr.value:>5.5g}"
                if attr.unit is None:
                    attr.unit = "None"
                if attr.description is None:
                    attr.description = ""

                param_table += f"{attr.name:>12} \t {value:>12} \t" +\
                    f"{attr.unit:>7} \t {attr.description} \n"

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
                    if isinstance(kwargs[arg], numbers.Number):
                        attr.value = kwargs[arg]     
                    elif isinstance(kwargs[arg], Parameter):
                        attr.value = kwargs[arg].value
                        attr.axis = kwargs[arg].axis

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
    

    def _to_dict(self):
        to_return = {"version": __version__, "type": "Pulse", "subclass": str(type(self))}

        for key, var in vars(self).items():
            if isinstance(var, Parameter):
                to_return[key] = var._to_dict()
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
                else:
                    return json.JSONEncoder.default(self, obj)
        
        return json.dumps(self._to_dict(), cls=autoEPREncoder, indent=4)

    def save(self, filename):
        """Save the Pulse to a JSON file.

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
        >>> obj = Pulse()
        >>> obj.save("my_pulse.json")
        """
                
        with open(filename, "w") as f:
           f.write(self._to_json())

    @classmethod
    def _from_dict(cls, dct):
        new_pulse = cls(tp=Parameter._from_dict(dct["tp"]), name=dct["name"])
        for key, var in dct.items(): 
            if isinstance(var, dict) and ("type" in var):
                setattr(new_pulse, key, Parameter._from_dict(var))
            elif key == "type":
                continue
            elif key == "version":
                continue
            else:
                setattr(new_pulse, key, var)
            
        return new_pulse

    @classmethod
    def _from_json(cls, JSONstring):
        dct = json.loads(JSONstring, object_hook=autoEPRDecoder)
        return cls._from_dict(dct)
    
    @classmethod
    def load(cls, filename):
        """Load a Pulse object from a JSON file.

        Parameters
        ----------
        filename : str
            Path to the JSON file.

        Returns
        -------
        obj : Pulse
            The Pulse loaded from the JSON file.

        Raises
        ------
        FileNotFoundError
            If the file does not exist.

        Example
        -------
        >>> obj = Pulse.load("my_pulse.json")
        """
        with open(filename, "r") as f:
           file_buffer = f.read()
        return cls._from_json(file_buffer)

# =============================================================================

#                                Subclasses

# =============================================================================

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
        Pulse.__init__(self, t=t, tp=tp, **kwargs)
        self.scale = None
        self.freq = Parameter("freq", freq, "GHz", "The detection frequency, if supported")
        self.pcyc = None
        pass


# =============================================================================


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
        self.pcyc = None
        self.scale = None

# =============================================================================


class RectPulse(Pulse):
    """
    Represents a rectangular monochromatic pulse.
    """

    def __init__(
            self, tp=16, freq=0, t=None, flipangle=None, **kwargs) -> None:
        """

        Parameters
        ----------
        tp : flaot, optional
            Pulse Length in ns, by default 16
        freq : float,optional
            Frequency in MHz, by default 0
        t : float, optional
            Time position in ns, by default None
        flipangle : _type_, optional
            The flip angle in radians, by default None
        """
        Pulse.__init__(
            self, tp=tp, t=t, flipangle=flipangle,name='RectPulse', **kwargs)
        self.freq = Parameter("freq", freq, "GHz", "Frequency of the Pulse")
        self._buildFMAM(self.func)
        pass

    def func(self, ax):
        nx = ax.shape[0]
        AM = np.ones(nx)
        FM = np.zeros(nx) + self.freq.value
        return AM, FM

class GaussianPulse(Pulse):
    """
    Represents a Gaussian monochromatic pulse.
    """

    def __init__(self, *, tp=32,FWHM=16, freq=0, **kwargs) -> None:
        """    Represents a Gaussian monochromatic pulse.

        Parameters
        ----------
        tp : float
            Pulse length in ns, by default 128
        FWHM : float,
            The full width at half maximum of the pulse
        freq : float, optional
            The frequency of the pulse, by default 0
        """
        Pulse.__init__(self, tp=tp,name='GaussianPulse', **kwargs)
        self.freq = Parameter("Freq", freq, "GHz", "Frequency of the Pulse")
        self.FWHM = Parameter("FWHM", FWHM, "ns", "Full Width at Half Maximum of the Pulse")
        self._buildFMAM(self.func)
        pass

    def func(self, ax):
        nx = ax.shape[0]
        sigma = self.FWHM.value / (2 * np.sqrt(2 * np.log(2)))
        AM = np.exp(-ax**2 / (2 * sigma**2))
        FM = np.zeros(nx) + self.freq.value
        return AM, FM
# =============================================================================
class FrequencySweptPulse(Pulse):
    """
    A general parent class for Frequency Swept Pulses.
    """
    
    def __init__(self, *, tp, t=None, scale=None, flipangle=None, pcyc=[0], name=None, **kwargs) -> None:
        
        super().__init__(tp=tp, t=t, scale=scale, flipangle=flipangle, pcyc=pcyc, name=name, **kwargs)

        # Extract frequency infomation, can be either specified with an initial and final frequency or a bandwidth combined with one of the two or the central frequency

        if "BW" in kwargs:
            BW = kwargs["BW"]
            if "init_freq" in kwargs:
                self.init_freq = Parameter("init_freq", kwargs["init_freq"], "GHz", "Initial frequency of pulse")
                self.final_freq = Parameter("final_freq", self.init_freq.value + BW, "GHz", "Final frequency of pulse")
            elif "final_freq" in kwargs:
                self.final_freq = Parameter("final_freq", kwargs["final_freq"], "GHz", "Final frequency of pulse")
                self.init_freq = Parameter("init_freq", self.final_freq.value - BW, "GHz", "Initial frequency of pulse")
            else:
                raise ValueError("Bandwidth must be combined with either an initial or final frequency")
        elif ("init_freq" in kwargs) & ("final_freq" in kwargs):
            self.init_freq = Parameter("init_freq", kwargs["init_freq"], "GHz", "Initial frequency of pulse")
            self.final_freq = Parameter("final_freq", kwargs["final_freq"], "GHz", "Final frequency of pulse")
        else:
            raise ValueError("Frequency information must be provided")
        
        

    @property
    def Qcrit(self):
        """ The critical Q factor for the pulse."""
        Qcrit = (2/np.pi)*np.log(2/(1+np.cos(self.flipangle.value)))
        Qcrit = np.min([Qcrit,5])
        return Parameter("Qcrit", Qcrit, None, "Critical Q factor for the pulse")
    
    @property
    def sweeprate(self):
        """ The sweep rate of the pulse in GHz/ns"""
        raise NotImplementedError("This property must be implemented in the subclass")

    @property
    def amp_factor(self):
        """ The B1 amplitude factor (nutation frequency) for the pulse in GHz"""
        amp_factor_value=  np.sqrt(2*np.pi*self.Qcrit.value*self.sweeprate.value)/(2*np.pi)
        return Parameter("amp_factor", amp_factor_value, "GHz", "Amplitude factor for the pulse")
    

class HSPulse(FrequencySweptPulse):
    """
    Represents a hyperboilc secant frequency-swept pulse.
    """
    def __init__(self, *, tp=128, order1=1, order2=6, beta=20, **kwargs) -> None:
        FrequencySweptPulse.__init__(self, tp=tp, name='HSPulse', **kwargs)
        self.order1 = Parameter(
            "order1", order1, None, "Order 1 of the HS Pulse")
        self.order2 = Parameter(
            "order2", order2, None, "Order 2 of the HS Pulse")
        self.beta = Parameter("beta", beta, None, "Beta of the HS Pulse")

        self._buildFMAM(self.func)
        pass


    def func(self, ax):
        beta = self.beta.value
        order1 = self.order1.value
        order2 = self.order2.value
        tp = ax.max() - ax.min()
        tcent = tp / 2
        ti = ax
        nx = ax.shape[0]
        # beta_exp1 = np.log(beta*0.5**(1-order1)) / np.log(beta)
        # beta_exp2 = np.log(beta*0.5**(1-order2)) / np.log(beta)
        # cut = round(nx/2)
        # AM = np.zeros(nx)
        # AM[0:cut] = 1/np.cosh(
        #     beta**beta_exp1 * (ax[0:cut]/tp)**order1)
        # AM[cut:-1] = 1/np.cosh(
        #     beta**beta_exp2 * (ax[cut:-1]/tp)**order2)

        # FM = BW * cumulative_trapezoid(AM**2,ax,initial=0) /\
        #      np.trapz(AM**2,ax) + self.init_freq.value
        sech = lambda x: 1/np.cosh(x) 
        cut = round(nx/2)
        AM = np.zeros_like(ti)
        AM[:cut] = sech(beta*0.5*(2*ti[:cut]/tp)**order1)
        AM[cut:] = sech(beta*0.5*(2*ti[cut:]/tp)**order2)


        BWinf = (self.bandwidth.value) / np.tanh(beta/2)

        freq = (BWinf/2) * np.tanh((beta/tp*ti)) + np.mean([self.init_freq.value,self.final_freq.value])
        # phase = 2*np.pi*(BWinf/2) * (tp/beta) * np.log(np.cosh((beta/tp)*ti))

        # total_phase = phase * 2* np.pi * np.mean([self.init_freq.value,self.final_freq.value])
        
        return AM, freq
    
    @property
    def sweeprate(self):
        """ The sweep rate of the pulse in GHz/ns"""
        sweeprate_value = self.beta.value * self.bandwidth.value / (2*self.tp.value)
        return Parameter("sweeprate", sweeprate_value, "GHz/ns", "Sweep rate of the pulse")

# =============================================================================

class ChirpPulse(FrequencySweptPulse):
    """
    Represents a linear frequency-swept pulse.
    """

    def __init__(self, *, tp=128, **kwargs) -> None:
        FrequencySweptPulse.__init__(self, tp=tp,name='ChirpPulse', **kwargs)

        self._buildFMAM(self.func)
        pass

    def func(self, ax):
        nx = ax.shape[0]
        AM = np.ones(nx)

        FM = np.linspace(
            self.init_freq.value, self.final_freq.value, nx)

        return AM, FM
    
    @property
    def sweeprate(self):
        """ The sweep rate of the pulse in GHz/ns"""
        sweeprate_value = self.bandwidth.value / self.tp.value
        return Parameter("sweeprate", sweeprate_value, "GHz/ns", "Sweep rate of the pulse")
# =============================================================================

class SincPulse(Pulse):


    def __init__(self, *, tp=128, freq=0, order=6, window=None, **kwargs) -> None:
        """    Represents a sinc shaped monochromatic pulse.


        Parameters
        ----------
        tp : int
            Pulse length in ns, by default 128
        freq : int, optional
            The frequency of the pulse, by default 0
        order : int, optional
            The order of this sinc function, by default 6
        window : _type_, optional
            The window function, by default None
        """

        Pulse.__init__(self, tp=tp,name='SincPulse', **kwargs)
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
    
# class GaussianPulse(Pulse):
#     """
#     Represents a Gaussian monochromatic pulse.
#     """

#     def __init__(self, *, tp=128, freq=0, **kwargs) -> None:
#         """    Represents a Gaussian monochromatic pulse.


#         Parameters
#         ----------
#         tp : int
#             Pulse length in ns, by default 128
#         freq : int, optional
#             The frequency of the pulse, by default 0
#         """

#         Pulse.__init__(self, tp=tp,name='GaussianPulse', **kwargs)
#         self.freq = Parameter("Freq", freq, "GHz", "Frequency of the Pulse")

#         self._buildFMAM(self.func)
#         pass


#     def func(self, ax):
#         nx = ax.shape[0]
#         FM = np.zeros(nx) + self.freq.value
#         AM = np.exp(-ax**2/(((self.tp.value/2)**2)/(-1*np.log(0.001))))

#         return AM, FM
    
# =============================================================================

def build_default_pulses(AWG=True,tp=12):

    exc_pulse = RectPulse(  
                    tp=tp, freq=0, flipangle=np.pi/2, scale=0
                )
    
    ref_pulse = RectPulse(  
                    tp=tp, freq=0, flipangle=np.pi, scale=0
                        )
    
    if AWG:
        pump_pulse = HSPulse(  
                    tp=120, init_freq=-0.25, final_freq=-0.03, flipangle=np.pi, scale=0,
                        order1=6, order2=1, beta=10
                    )
        det_event = Detection(tp=512, freq=0)
    else:
        pump_pulse = RectPulse(
                    tp=tp, freq=-0.07, flipangle=np.pi, scale=0)
        
        # det_event = Detection(tp=exc_pulse.tp * 2, freq=0)
        det_event = Detection(tp=128, freq=0)
        
    pulses = {'exc_pulse':exc_pulse, 'ref_pulse':ref_pulse, 'pump_pulse':pump_pulse, 'det_event':det_event}

    return pulses