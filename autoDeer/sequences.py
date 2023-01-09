from autodeer.hardware import Sequence, RectPulse, Pulse, ChirpPulse, \
    Detection, Delay, Parameter
import numpy as np

class DEERSequence(Sequence):

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
        super().__init__(name, B, LO, reptime, averages, shots, **kwargs)

    def four_pulse(self, tp):

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
            t=self.tau1 + self.deadtime, tp=tp, freq=0, flipangle=np.pi
        ))
        self.addPulse(RectPulse( # Ref 2 pulse
            t=2*self.tau1 + self.tau2, tp=tp, freq=0, flipangle=np.pi
        ))
        self.addPulse(Detection(t=2*(self.tau1+self.tau2), tp=512))

        self.addPulsesProg(
            2,
            ["t"],
            0,
            axis)


        pass
    def five_pulse(self, tp):

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
            3,
            ["t"],
            0,
            axis)
        pass

    def select_pcyc(self, option: str):
        """Choose which DEER phase you would like.

        .. |xp| replace:: x\ :sub:`p`\

        +---------------+-------------+----------------+--------+---------------------------+-----------------------------+------------+
        | Phase cycle   | Short Code  | DEER_sequence  | Steps  | Pulse Phase Cycle         | Remaining Echoes            | Reference  |
        +===============+=============+================+========+===========================+=============================+============+
        | (x)x|xp|x     | DC          | ALL            | 2      | [+(+x) -(-x)]             | PE12rp,SE(PE12)p3,PE12rpr3  |            |
        +---------------+-------------+----------------+--------+---------------------------+-----------------------------+------------+
        | x[x][|xp|]x   | 16step_4p   | 4 pulse        | 16     | [+(+x) -(+y)+(-x) -(-y)]  | -                           | [1]        |
        +---------------+-------------+----------------+--------+                           +-----------------------------+------------+
        | xxp[x][|xp|]x | 16step_5p   | 5 pulse        | 16     | [+(+x) +(+y)+(-x) +(-y)]  | PEp02r3,b PE1p0r2r3b        | [1]        |
        +---------------+-------------+----------------+--------+---------------------------+-----------------------------+------------+


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

    def __init__(self, *, B, LO, reptime, averages, shots, **kwargs) -> None:
        
        name = "Hahn Echo"

        super().__init__(name, B, LO, reptime, averages, shots, **kwargs)

        tau = 500
        tp = 12

        name = "Hahn Echo"
        self.addPulse(RectPulse(  # Exc pulse
            t=0, tp=tp, freq=0, flipangle=np.pi/2, pcyc=([0, np.pi], [1 -1])
        ))
        self.addPulse(RectPulse( # Pump 1 pulse
            t=tau, tp=tp, freq=0, flipangle=np.pi
        ))

        self.addPulse(Detection(t=2*tau, tp=self.det_window))

class FieldSweepSequence(HahnEchoSequence):

    def __init__(self, *, B, LO, Bwidth, reptime, averages, shots, **kwargs) -> None:
        
        super().__init__(B, LO, reptime, averages, shots, **kwargs)
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

    def __init__(self, *, name, B, LO, reptime, averages, shots, **kwargs) -> None:
        
        name = "Carr-Purcell"
        super().__init__(name, B, LO, reptime, averages, shots, **kwargs)