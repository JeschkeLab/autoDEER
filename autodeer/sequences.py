import numpy as np
from autodeer import __version__
from pyepr import Sequence, Parameter, Detection, RectPulse, val_in_us, add_phaseshift, val_in_ns
import deerlab as dl

class DEERSequence(Sequence):
    
    """
    Represents a DEER/PELDOR sequence. 
    """

    def __init__(
        self, *, tau1, tau2, tau3 = None, tau4=None, dt, B, freq, reptime,
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
        freq : int or float
            The freq frequency in GHz
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
        name = "DEERSequence"
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
        self.relaxation = False
        self.add_ESEEM_avg(None)

        super().__init__(
            name=name, B=B, freq=freq, reptime=reptime, averages=averages,
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
            self.relaxation = True
        else:
            dim = np.floor((self.tau2.value)/self.dt)
            if dim < 2:
                print(f"dim: {dim}, tau2: {self.tau2.value}, dt: {self.dt}")
                raise ValueError("Time axis is too short, increase tau2")
            self.t = Parameter(name="t", value = self.tau1.value - self.deadtime, step=self.dt,
                       dim=dim, unit="ns", description="The time axis", virtual=True)
            self.relaxation = False

        
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
            self.relaxation = True

        else:
            if self.deadtime > self.tau3.value:
                raise ValueError("Deadtime must be greater than tau3, otherwise pulses crash")
            
            dim = np.floor((self.tau2.value + self.tau1.value - self.tau3.value)/self.dt)
            self.t = Parameter(
                name="t", value= self.tau3.value - self.deadtime, step=self.dt,
                dim=dim, unit="ns", description="The time axis", virtual=True)
            self.relaxation = False

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
            self.relaxation = True

        else:
            dim = np.floor((self.tau2.value + self.tau1.value - 2*self.deadtime)/self.dt)
            self.t = Parameter(
                name="t", value = + self.deadtime, step=self.dt,
                dim=dim, unit="ns", description="The time axis", virtual=True)
            self.relaxation = False

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
            self.relaxation = True

        else:
            
            dim = np.floor((n*self.step.value)/self.dt)

            self.t = Parameter(
                name="t", value= - self.deadtime, step=self.dt,
                dim=dim, unit="ns", description="The time axis", virtual=True)
            self.relaxation = False
                

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

        self.pcyc_name = option
        
        self._buildPhaseCycle()


    def simulate(self):
        """"
        Simulate the DEER sequence using a DeerLab dipolar model for a gaussian of mean 4.5 nm and std 1.0 nm, with a spin concentration of 50 uM.
        
        Returns
        -------
        t : np.ndarray
            The time axis in us
        Vsim : np.ndarray
            The simulated DEER trace
        """
        if self.relaxation:

            T_CP=2.5e3
            e=1.8
            xaxis = val_in_ns(self.tau2)
            func = lambda x, a, tau, e: a*np.exp(-(x/tau)**e)
            data = func(xaxis,1,T_CP,e)
            data = add_phaseshift(data, 0.05)
            return xaxis, data
        else:

            if self.name == "4pDEER":
                exp_type = "4pDEER"
                tau1 = val_in_us(self.tau1)
                tau2 = val_in_us(self.tau2)
                t = val_in_us(self.t)
            elif self.name == "5pDEER":
                exp_type = "5pDEER"
                tau1 = val_in_us(self.tau1)
                tau2 = val_in_us(self.tau2)
                tau3 = val_in_us(self.tau3)
                t = val_in_us(self.t)
            elif self.name == "3pDEER":
                exp_type = "3pDEER"
                tau1 = val_in_us(self.tau1)
                t = val_in_us(self.t)
            elif self.name == "nDEER-CP":
                exp_type = "4pDEER"
                tau1 = val_in_us(self.tau1)
                tau2 = val_in_us(self.tau2)
                t = val_in_us(self.t)

            if exp_type == "4pDEER":
                experimentInfo = dl.ex_4pdeer(tau1=tau1,tau2=tau2,pathways=[1,2,3])
                reftimes = dict(zip(["reftime1","reftime2","reftime3"],experimentInfo.reftimes(tau1,tau2)))
                mod_depths = {"lam1":0.4, "lam2":0.1, "lam3":0.2}
            elif exp_type == "5pDEER":
                experimentInfo = dl.ex_fwd5pdeer(tau1=tau1,tau2=tau2,tau3=tau3,pathways=[1,2,3,4,5])
                reftimes = dict(zip(["reftime1","reftime2","reftime3","reftime4","reftime5"],experimentInfo.reftimes(tau1,tau2,tau3)))
                mod_depths = {"lam1":0.4, "lam2":0.00, "lam3":0.0, "lam4":0.00, "lam5":0.1}

            elif exp_type == "3pDEER":
                experimentInfo = dl.ex_3pdeer(tau=tau1,pathways=[1,2])
                reftimes = dict(zip(["reftime1","reftime2"],experimentInfo.reftimes(tau1,)))
                mod_depths = {"lam1":0.6, "lam2":0.1}

            r = np.linspace(0.5,10,100)
            rmean = 4.5
            rstd = 1.0
            conc = 50

            Vmodel = dl.dipolarmodel(t,r,Pmodel=dl.dd_gauss, experiment=experimentInfo)
            Vsim = Vmodel(mean=rmean, std=rstd, conc=conc, scale=1, **reftimes, **mod_depths)
            # Add phase shift
            Vsim = add_phaseshift(Vsim, 0.05)
            return t, Vsim

# =============================================================================

class RefocusedEcho2DSequence(Sequence):

    """
    Represents a 2D Refocused-echo Sequence. 
    """
    def __init__(self, *, B, freq, reptime, averages, shots,
            tau, dim=100, **kwargs) -> None:
        """Build a 2D Refocused-echo sequence using either 
        rectangular pulses or specified pulses.

        Parameters
        ----------
        B : int or float
            The B0 field, in Guass
        freq : int or float
            The freq frequency in GHz
        reptime : _type_
            The shot repetition time in us
        averages : int
            The number of scans.
        shots : int
            The number of shots per point
        tau : int
            The maximum total sequence length in us
        dim: int
            The number of points in both the X and Y axis
    

        Optional Parameters
        -------------------
        pi2_pulse : Pulse
            An autoEPR Pulse object describing the excitation pi/2 pulse. If
            not specified a RectPulse will be created instead. 
        pi_pulse : Pulse
            An autoEPR Pulse object describing the refocusing pi pulses. If
            not specified a RectPulse will be created instead. 
        """

        name = "RefocusedEcho2D"
        super().__init__(
            name=name, B=B, freq=freq, reptime=reptime, averages=averages,
            shots=shots, **kwargs)
        
        tau_init = 400
        dt = (tau * 1e3) / dim

        self.tau1 = Parameter(name="tau1", value=tau_init, dim=dim, step=dt, unit="ns",
            description="1st interpulse delay",virtual=True)
        self.tau2 = Parameter(name="tau2", value=tau_init, dim=dim, step=dt, unit="ns",
            description="2nd interpulse delay",virtual=True)


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

    def simulate(self,sigma = 0.8):
        """
        Simulates the Refocused-echo 2D sequence using this function.
        .. math::
            V(tau1, tau2) = e^{-(tau1^2 + tau2^2 - tau1*tau2)/(2*sigma^2)}
        
        Parameters
        ----------
        sigma : float
            The decay rate, by default 0.8 us
        
        Returns
        -------
        tau1, tau2 : list[np.ndarray]
            The two time axis in us
        data : np.ndarray
            The simulated 2D data
        """
        func = lambda x, y: np.exp(-((x**2 + y**2 - 1*x*y) / (2*sigma**2)))

        xaxis = val_in_us(self.tau1)
        yaxis = val_in_us(self.tau2)

        X, Y = np.meshgrid(xaxis, yaxis)
        data = func(X, Y)
        data = add_phaseshift(data, 0.05)
        return [xaxis, yaxis], data

