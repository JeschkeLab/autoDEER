from autodeer.openepr import Sequence, RectPulse, Pulse, ChirpPulse, \
    Detection, Delay, Parameter
import numpy as np

class DEERSequence(Sequence):
    """
    Represents a DEER/PELDOR sequence. 
    """

    def __init__(
        self, *, tau1, tau2, tau3 = None, dt, B, LO, reptime, averages, shots,
        **kwargs) -> None:
        """Build a DEER sequence using rectangular pulses

        Parameters
        ----------
        tau1 : int or float
            The first interpulse delay in ns
        tau2 : int or float
            The second interpulse delay in ns
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
            the 1st refocusing pulse in ns, by default None. If the value is 
            None then a 4-pulse sequence is created, otherwise a 5-pulse. 
        
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
        self.tau1 = tau1
        self.tau2 = tau2
        self.tau3 = tau3
        self.dt = dt
        self.deadtime = 200
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

    def four_pulse(self, tp=16):
        """
        Build a four pulse DEER sequence.

        Parameters
        ----------
        tp : float
            Length of default RectPulse in ns, by default 16ns.
        """

        axis = np.arange(
            2*self.tau1-self.deadtime, self.tau2 + 2*self.tau1 - self.deadtime,
            self.dt)

        if hasattr(self,"exc_pulse"): # Exc pulse
            self.addPulse(self.exc_pulse.copy(t=0))
        else:
            self.addPulse(RectPulse(  
                t=0, tp=tp, freq=0, flipangle=np.pi/2
            ))
        
        if hasattr(self,"ref_pulse"): # Ref 1 pulse
            self.addPulse(self.ref_pulse.copy(t=self.tau1))
        else:
            self.addPulse(RectPulse( 
                t=self.tau1, tp=tp, freq=0, flipangle=np.pi
            ))

        if hasattr(self,"pump_pulse"): # Pump 1 pulse
            self.addPulse(self.pump_pulse.copy(t=2*self.tau1 - self.deadtime))
        else:
            self.addPulse(RectPulse(
                t=2*self.tau1 - self.deadtime, tp=tp, freq=0, flipangle=np.pi
            ))

        if hasattr(self,"ref_pulse"): # Ref 2 pulse
            self.addPulse(self.ref_pulse.copy(t=2*self.tau1 + self.tau2))
        else:
            self.addPulse(RectPulse( 
                t=2*self.tau1 + self.tau2, tp=tp, freq=0, flipangle=np.pi
            ))

        if hasattr(self, "det_event"):
            self.addPulse(self.det_event.copy(t=2*(self.tau1+self.tau2)))
        else:
            self.addPulse(Detection(t=2*(self.tau1+self.tau2), tp=512))

        self.addPulsesProg(
            [2],
            ["t"],
            0,
            axis)


        pass
    def five_pulse(self, tp=16):
        """
        Build a five pulse DEER sequence.

        Parameters
        ----------
        tp : float
            Length of default RectPulse in ns, by default 16ns.
        """
        axis = np.arange(
            self.tau1+self.deadtime, self.tau2 + 2*self.tau1 - self.deadtime,
             self.dt)


        if hasattr(self,"exc_pulse"): # Exc pulse
            self.addPulse(self.exc_pulse.copy(t=0))
        else:
            self.addPulse(RectPulse(  
                t=0, tp=tp, freq=0, flipangle=np.pi/2
            ))

        if hasattr(self,"pump_pulse"): # Pump 1 pulse
            self.addPulse(self.pump_pulse.copy(t=self.tau1 - self.tau3))
        else:
            self.addPulse(RectPulse(
                t=self.tau1 - self.tau3, tp=tp, freq=0, flipangle=np.pi
            ))

        if hasattr(self,"ref_pulse"): # Ref 1 pulse
            self.addPulse(self.ref_pulse.copy(t=self.tau1))
        else:
            self.addPulse(RectPulse( 
                t=self.tau1, tp=tp, freq=0, flipangle=np.pi
            ))

        if hasattr(self,"pump_pulse"): # Pump 2 pulse
            self.addPulse(self.pump_pulse.copy(t=self.tau1 + self.deadtime))
        else:
            self.addPulse(RectPulse(
                t=self.tau1 + self.deadtime, tp=tp, freq=0, flipangle=np.pi
            ))

        if hasattr(self,"ref_pulse"): # Ref 2 pulse
            self.addPulse(self.ref_pulse.copy(t=2*self.tau1 + self.tau2))
        else:
            self.addPulse(RectPulse( 
                t=2*self.tau1 + self.tau2, tp=tp, freq=0, flipangle=np.pi
            ))

        if hasattr(self, "det_event"):
            self.addPulse(self.det_event.copy(t=2*(self.tau1+self.tau2)))
        else:
            self.addPulse(Detection(t=2*(self.tau1+self.tau2), tp=512))

        self.addPulsesProg(
            [3],
            ["t"],
            0,
            axis)
        pass

    def seven_pulse(self, tp):
        """
        Build a seven pulse DEER sequence.

        Parameters
        ----------
        tp : float
            Length of default RectPulse in ns, by default 16ns.
        """
        self.name = "7p-DEER"

        axis = np.arange(
            self.tau1+self.deadtime, self.tau2 + 2*self.tau1 - self.deadtime,
             self.dt)
        self.addPulse(RectPulse(  # Exc pulse
            t=0, tp=tp, freq=0, flipangle=np.pi/2
        ))
        self.addPulse(RectPulse( # Pump 1 pulse
            t=self.tau1 - self.tau3, tp=tp, freq=0, flipangle=np.pi
        ))

        self.addPulse(RectPulse( # Ref 1 pulse
            t=self.tau1, tp=tp, freq=0, flipangle=np.pi
        ))
        self.addPulse(RectPulse( # Pump 2 pulse
            t=self.tau1 + self.deadtime, tp=tp, freq=0, flipangle=np.pi
        ))
        self.addPulse(RectPulse( # Ref 2 pulse
            t=2*self.tau1 + self.tau2, tp=tp, freq=0, flipangle=np.pi
        ))
        self.addPulse(RectPulse( # Pump 3 pulse
            t=self.tau1 + self.deadtime, tp=tp, freq=0, flipangle=np.pi
        ))
        self.addPulse(RectPulse( # Ref 3 pulse
            t=2*self.tau1 + self.tau2, tp=tp, freq=0, flipangle=np.pi
        ))
        self.addPulse(Detection(t=2*(self.tau1+self.tau2), tp=512))

        self.addPulsesProg(
            [3],
            ["t"],
            0,
            axis)
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
        
        self._buildPhaseCycle()


    def simulate(self):
        t = np.arange(0, self.tau1+self.tau2, self.dt)
        r = np.arange(1.5, 8)
        pass


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

        tau = 500
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
            self.addPulse(self.pi_pulse.copy(t=tau))
        else:
            self.addPulse(RectPulse( # Pump 1 pulse
                t=tau, tp=tp, freq=0, flipangle=np.pi
            ))

        if hasattr(self, "det_event"):
            self.addPulse(self.det_event.copy(t=2*tau))
        else:
            self.addPulse(Detection(t=2*tau, tp=self.det_window.value))

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


        self.Bwidth = Parameter(
            "Bwidth", Bwidth, "Gauss", "Field sweep width"
        )

        self.addPulseProg(
            pulse_id=None,
            variable='B',
            axis_id=0,
            axis = np.linspace(B-Bwidth/2, B+Bwidth/2, Bwidth + 1)
        )
        

class CarrPurcellSequence(Sequence):
    """
    Represents a Carr-Purcell sequence. 
    """
    def __init__(self, *, B, LO, reptime, averages, shots,
            tau, n, **kwargs) -> None:
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
            The maximum total sequence length in ns
        n : int
            The number refocusing pulses

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
        self.tau = Parameter(name="tau", value=tau, unit="ns",
            description="Total sequence length")
        self.n = Parameter(name="n", value=n,
            description="The number of pi pulses", unit="None")

        if "pi_pulse" in kwargs:
            self.pi_pulse = kwargs["pi_pulse"]
        if "pi2_pulse" in kwargs:
            self.pi2_pulse = kwargs["pi2_pulse"]

        self._build_sequence()

    def _build_sequence(self):

        tau = self.tau.value
        n = self.n.value
        step = tau/(2*n)
        multipliers = [1]
        multipliers += [1+2*i for i in range(1,n)]
        multipliers += [2*n]
        deadtime = 300

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
        
        self.addPulse(Detection(t=tau, tp=512))
        
        axis = np.arange(deadtime,tau/(2*n),10)
        self.addPulsesProg(
            list(range(1,n+2)),
            ["t"]*(n+1),
            0,
            axis,
            multipliers=multipliers)



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

        self.addPulse(RectPulse(  # Hard pulse
            t=0, tp=4, freq=0, flipangle="Hard"
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
        nut_axis = np.arange(0,66,2,)
        self.addPulsesProg(
            [0],
            ["tp"],
            0,
            nut_axis)

        # Add frequency sweep
        width= 0.3
        axis = np.arange(self.LO.value-width,self.LO.value+width+0.02,0.02)
        self.addPulsesProg(
            [None, None],
            ["LO", "B"],
            1,
            axis,
            multipliers=[1,1/self.gyro])


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
