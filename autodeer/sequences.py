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
            The first interpulse delay
        tau2 : int or float
            The second interpulse delay
        dt : int or float
            The time step for DEER measurment
        B : int or float
            The B0 field, in Guass
        LO : int or float
            The LO frequency
        reptime : _type_
            _description_
        averages : int
            The number of scans.
        shots : int
            The number of shots per point
        tau3 : int or float, optional
            The delay between the first static pump pulse in 5-pulse DEER and 
            the 1st refocusing pulse, by default None. If the value is None
            then a 4-pulse sequence is created, otherwise a 5-pulse. 
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

    def four_pulse(self, tp):
        """
        Build a four pulse DEER sequence.

        Parameters
        ----------
        tp : float
            Step size in ns
        """

        axis = np.arange(
            2*self.tau1-self.deadtime, self.tau2 + 2*self.tau1 - self.deadtime,
            self.dt)

        self.addPulse(RectPulse(  # Exc pulse
            t=0, tp=tp, freq=0, flipangle=np.pi/2
        ))
        self.addPulse(RectPulse( # Ref 1 pulse
            t=self.tau1, tp=tp, freq=0, flipangle=np.pi
        ))
        self.addPulse(RectPulse( # Pump 1 pulse
            t=2*self.tau1 - self.deadtime, tp=tp, freq=0, flipangle=np.pi
        ))
        self.addPulse(RectPulse( # Ref 2 pulse
            t=2*self.tau1 + self.tau2, tp=tp, freq=0, flipangle=np.pi
        ))
        self.addPulse(Detection(t=2*(self.tau1+self.tau2), tp=512))

        self.addPulsesProg(
            [2],
            ["t"],
            0,
            axis)


        pass
    def five_pulse(self, tp):
        """
        Build a five pulse DEER sequence.

        Parameters
        ----------
        tp : float
            Step size in ns
        """
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
            Step size in ns
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
            self.pulses[1]._addPhaseCycle(
                [0, np.pi/2, np.pi, -np.pi/2], [1, -1, 1, -1])
            self.pulses[2]._addPhaseCycle(
                [0, np.pi/2, np.pi, -np.pi/2], [1, 1, 1, 1])
            
        elif option == "16step_5p":
            self.pulses[2]._addPhaseCycle(
                [0, np.pi/2, np.pi, -np.pi/2], [1, -1, 1, -1])
            self.pulses[3]._addPhaseCycle(
                [0, np.pi/2, np.pi, -np.pi/2], [1, 1, 1, 1])

    def simulate(self):
        t = np.arange(0, self.tau1+self.tau2, self.dt)
        r = np.arange(1.5, 8)
        pass


class HahnEchoSequence(Sequence):
    """
    Represents a Hahn-Echo sequence. 
    """
    def __init__(self, *, B, LO, reptime, averages, shots, **kwargs) -> None:
        
        name = "Hahn Echo"

        super().__init__(
            name=name, B=B, LO=LO, reptime=reptime, averages=averages,
            shots=shots, **kwargs)

        tau = 500
        tp = 12

        name = "Hahn Echo"
        self.addPulse(RectPulse(  # Exc pulse
            t=0, tp=tp, freq=0, flipangle=np.pi/2, pcyc={"phases":[0, np.pi], "dets":[1, -1]}
        ))
        self.addPulse(RectPulse( # Pump 1 pulse
            t=tau, tp=tp, freq=0, flipangle=np.pi
        ))

        self.addPulse(Detection(t=2*tau, tp=self.det_window.value))

class FieldSweepSequence(HahnEchoSequence):
    """
    Represents a Field Sweep (EDFS) sequence. 
    """
    def __init__(self, *, B, LO, Bwidth, reptime, averages, shots, **kwargs) -> None:
        
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
    def __init__(self, *, B, LO, reptime, averages, shots, **kwargs) -> None:
        
        name = "Carr-Purcell"
        super().__init__(
            name=name, B=B, LO=LO, reptime=reptime, averages=averages,
            shots=shots, **kwargs)

class ResonatorProfileSequence(Sequence):
    """
    Builds nutation based Resonator Profile sequence. 
    """

    def __init__(self,*,B,LO,reptime,averages,shots,**kwargs) -> None:

        name = "Resonator-Profile"
        super().__init__(
            name=name, B=B, LO=LO, reptime=reptime, averages=averages,
            shots=shots, **kwargs)
        self.gyro = LO/B

        self._build_sequence()

    def _build_sequence(self):

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
                [None],
                ["LO"],
                1,
                axis)
   
        self.addPulsesProg(
                [None],
                ["B"],
                1,
                axis/self.gyro)


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
