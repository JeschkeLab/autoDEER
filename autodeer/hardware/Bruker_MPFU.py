from autodeer.classes import Interface, Parameter
from autodeer.dataset import Dataset
from autodeer.pulses import Delay, Detection, RectPulse
from autodeer.hardware.XeprAPI_link import XeprAPILink
from autodeer.hardware.Bruker_tools import PulseSpel, run_general,PSPhaseCycle, write_pulsespel_file
from autodeer.sequences import Sequence, HahnEchoSequence
from autodeer.utils import save_file

import tempfile
import time
from scipy.optimize import minimize_scalar, curve_fit
import numpy as np

# =============================================================================


class BrukerMPFU(Interface):
    """
    Represents the interface for connecting to MPFU based Bruker ELEXSYS-II 
    Spectrometers.
    """
    def __init__(self, config_file:str, d0=600) -> None:
        """An interface for connecting to MPFU based Bruker ELEXSYS-II 
        Spectrometers.

        Getting Started
        ------------------
        Before a connection can be made an appropriate configuration file first
        needs to be written. 

        1. Open Xepr
        2. Processing -> XeprAPI -> Enable XeprAPI
        3. `BrukerAWG.connect()`

        Parameters
        ----------
        config_file : str
            _description_
        d0 : int, optional
            _description_, by default 600
        """

        self.api = XeprAPILink(config_file)
        self.spec_config = self.api.config["Spectrometer"]
        self.bridge_config = self.api.spec_config["Bridge"]

        self.MPFU = self.bridge_config["MPFU Channels"]
        
        self.temp_dir = tempfile.mkdtemp("autoDEER")

        self.d0 = d0
        
        super().__init__()


    def connect(self) -> None:

        self.api.connect()
        return super().connect()

    def acquire_dataset(self) -> Dataset:
        return self.api.acquire_dataset()
    
    def launch(self, sequence: Sequence, savename: str, start=True, tune=True):
                    
        channels = _MPFU_channels(sequence)
        # pcyc = PSPhaseCycle(sequence,self.MPFU)

        N_channels = len(channels)

        if N_channels > len(self.MPFU):
            raise RuntimeError(
                f"This sequence requires {N_channels} MPFU" "Channels." 
                "Only {len(self.MPFU)} are avaliable on this spectrometer.")
        
        if tune:
            MPFUtune(self,sequence, channels)

        # PSpel = PulseSpel(sequence, MPFU=self.MPFU)
        def_file, exp_file = write_pulsespel_file(sequence,False,self.MPFU)
        PSpel_file = self.temp_dir + "/autoDEER_PulseSpel"
        # PSpel.save(PSpel_file)
        save_file(PSpel_file +'.def',def_file)
        save_file(PSpel_file +'.exp',exp_file)
        self.api.set_field(sequence.B.value)
        self.api.set_freq(sequence.LO.value)
        if hasattr(sequence,"Bwidth"):
            self.api.set_sweep_width(sequence.Bwidth.value)
        
        run_general(self.api,
            ps_file= [PSpel_file],
            exp=("auto","auto"),
            settings={"ReplaceMode": False},
            variables={"d0": self.d0},
            run=False
        )
        if start:
            self.api.run_exp()
        pass
    
    def tune(self, sequence, B0, LO) -> None:
        channels = _MPFU_channels(sequence)
        
        for i,channel in enumerate(channels):
            ps_length = int(np.pi /(channel[0]*2))
            phase = channel[1]
            if phase == 0:
                echo = "R+"
            elif phase == np.pi/2:
                echo = "I+"
            elif phase == np.pi:
                echo = "R-"
            elif (phase == -np.pi/2) or (phase == 3*np.pi/2):
                echo = "I-"
            mpfu_tune = MPFUtune(self.api,B0=B0,LO=LO, echo="Hahn",ps_length=ps_length)
            mpfu_tune.tune({self.MPFU[i]: echo})
        
        pass

    def isrunning(self) -> bool:
        return self.api.is_exp_running()

    def terminate(self) -> None:
        return self.api.abort_exp()

# =============================================================================


def _MPFU_channels(sequence):
    """Idenitifies how many unique MPFU channels are needed for a sequence and
    applies the correct Channel infomation to each pulse.
    """
    channels = []

    for iD, pulse in enumerate(sequence.pulses):
        if type(pulse) is Delay:
            continue
        if type(pulse) is Detection:
            continue

        if pulse.tp.value == 0:
            flip_power = np.inf
        else:
            flip_power = pulse.flipangle.value / pulse.tp.value

        if (pulse.freq.value != 0):
            # This is an ELDOR pulse
            pulse.pcyc["Channels"] = "ELDOR"
            channels.append("ELDOR")
            continue

        if not "Channels" in pulse.pcyc:
            pulse.pcyc["Channels"] = []
        elif pulse.pcyc is None:
            pulse.pcyc = {}
            pulse.pcyc["Channels"] = []
        
        
        for phase in pulse.pcyc["Phases"]:
            power_phase = (flip_power, phase)
            if power_phase in channels:
                channel = channels.index(power_phase)
            else:
                channels.append(power_phase)
                channel = channels.index(power_phase)
            pulse.pcyc["Channels"].append(channel)
    
    return channels

# =============================================================================



def MPFUtune(interface, sequence, channels, echo='Hahn',tol: float = 0.1,
             bounds=[0, 100]):

    hardware_wait=5
    def tune_phase(
        channel: str, target: str, tol=0.1, maxiter=30) -> float:
        """Tunes the phase of a given channel to a given target using the
        standard scipy optimisation scripts. 

        Parameters
        ----------
        channel : str
            The chosen MPFU channel. Options: ['+<x>', '-<x>', '+<y>', '-<y>']
        target : str
            The target echo position, this can either be maximising (+) or
            minimising (-) either the real (R) or imaginary (I) of the echo. 
            Options: ['R+', 'R-', 'I+', 'I-']
        tol : float, optional
            The tolerance in phase parameter, by default 0.1
        maxiter : int, optional
            The maximum number of iterations in the optimisation, by default 30

        Returns
        -------
        float
            The optimal value of the phase parameter

        """

        channel_opts = ['+<x>', '-<x>', '+<y>', '-<y>']
        phase_opts = ['R+', 'R-', 'I+', 'I-']

        if channel not in channel_opts:
            raise ValueError(f'Channel must be one of: {channel_opts}')
    
        if target not in phase_opts:
            raise ValueError(f'Phase target must be one of: {phase_opts}')
        if channel == '+<x>':
            phase_channel = 'BrXPhase'
        elif channel == '-<x>':
            phase_channel = 'BrMinXPhase'
        elif channel == '+<y>':
            phase_channel = 'BrYPhase'
        elif channel == '-<y>':
            phase_channel = 'BrMinYPhase'

        if target == 'R+':
            test_fun = lambda x: -1 * np.real(x)
        elif target == 'R-':
            test_fun = lambda x: 1 * np.real(x)
        elif target == 'I+':
            test_fun = lambda x: -1 * np.imag(x)
        elif target == 'I-':
            test_fun = lambda x: 1 * np.imag(x)

        lb = 0.0
        ub = 100.0

        def objective(x, *args):
            # x = x[0]
            interface.api.hidden[phase_channel].value = x  # Set phase to value
            time.sleep(hardware_wait)
            interface.api.run_exp()
            while interface.api.is_exp_running():
                time.sleep(1)
            data = interface.api.acquire_scan()
            v = data.data

            val = test_fun(np.sum(v))

            print(f'Phase Setting = {x:.1f} \t Echo Amplitude = {-1*val:.2f}')

            return val

        output = minimize_scalar(
            objective, method='bounded', bounds=[lb, ub],
            options={'xatol': tol, 'maxiter': maxiter})
        result = output.x
        print(f"Optimal Phase Setting for {phase_channel} is: {result:.1f}")
        interface.api.hidden[phase_channel].value = result
        return result

    def tune_power(
            channel: str, tol=0.1, maxiter=30,
            bounds=[0, 100]) -> float:
            """Tunes the attenuator of a given channel to a given target using the
            standard scipy optimisation scripts. 

            Parameters
            ----------
            channel : str
                The chosen MPFU channel. Options: ['+<x>', '-<x>', '+<y>', '-<y>']
            tol : float, optional
                The tolerance in attenuator parameter, by default 0.1
            maxiter : int, optional
                The maximum number of iterations in the optimisation, by default 30

            Returns
            -------
            float
                The optimal value of the attenuator parameter

            """
            channel_opts = ['+<x>', '-<x>', '+<y>', '-<y>']
            if channel not in channel_opts:
                raise ValueError(f'Channel must be one of: {channel_opts}')
            
            if channel == '+<x>':
                atten_channel = 'BrXAmp'
            elif channel == '-<x>':
                atten_channel = 'BrMinXAmp'
            elif channel == '+<y>':
                atten_channel = 'BrYAmp'
            elif channel == '-<y>':
                atten_channel = 'BrMinYAmp'

            lb = bounds[0]
            ub = bounds[1]

            def objective(x, *args):
                interface.api.hidden[atten_channel].value = x  # Set phase to value
                time.sleep(hardware_wait)
                interface.api.run_exp()
                while interface.api.is_exp_running():
                    time.sleep(1)
                data = interface.api.acquire_scan()
                v = data.data

                val = -1 * np.sum(np.abs(v))

                print(f'Power Setting = {x:.1f} \t Echo Amplitude = {-1*val:.2f}')

                return val

            output = minimize_scalar(
                objective, method='bounded', bounds=[lb, ub],
                options={'xatol': tol, 'maxiter': maxiter})
            result = output.x
            print(f"Optimal Power Setting for {atten_channel} is: {result:.1f}")
            interface.api.hidden[atten_channel].value = result
            return result
    
    for i,channel in enumerate(channels):
        MPFU_chanel = interface.MPFU[i]
        if MPFU_chanel == '+<x>':
            phase_cycle = 'BrXPhase'
        elif MPFU_chanel == '-<x>':
            phase_cycle = 'BrMinXPhase'
        elif MPFU_chanel == '+<y>':
            phase_cycle = 'BrYPhase'
        elif MPFU_chanel == '-<y>':
            phase_cycle = 'BrMinYPhase'
        
        tau = Parameter("tau",tau,dim=4,step=0)

        exc_pulse = RectPulse(freq=0, tp = np.around(channel[0] * np.pi/2), scale=1, flipangle = np.pi/2)
        ref_pulse = RectPulse(freq=0, tp = np.around(channel[0] * np.pi), scale=1, flipangle = np.pi)
        
        seq = HahnEchoSequence(
            B=sequence.B,LO=sequence.LO,reptime=sequence.reptime,averages=1,
            shots=10, tau=tau, exc_pulse=exc_pulse, ref_pulse=ref_pulse)
        
        seq.evolution([tau])

        interface.launch(seq, savename="autoTUNE", start=False, tune=False)


        interface.api.set_PulseSpel_phase_cycling(phase_cycle)
        print(f"Tuning channel: {channel}")
        tune_power(channel, tol=tol, bounds=bounds)
        tune_phase(channel, channels[channel], tol=tol)




class MPFUtune1:
    """
    Tuning MPFU channels for optimal attenuation and phase
    """
    def __init__(
        self, interface, B, LO, echo="Hahn", ps_length=16, d0=680, srt=6e6,) -> None:
        """
        Parameters
        ----------
        api : autoEPR.Interface
            The spectrometer interface
        B0 : float
            The B0 magnetic field, in Gauss
        LO : float,
            The MW frequency, in GHz
        echo : str, optional
            The echo type. Options = ['Hahn","Refocused"], by default "Hahn"
        ps_length : int, optional
            The length of the pi/2 pulse, by default 16
        d0 : int, optional
            The approximate position of d0, this should be lower than ideal,
             by default 680
        """
        self.interface = interface
        self.api = self.interface.api
        self.MPFU = self.interface.MPFU
        self.hardware_wait = 5  # seconds 
        self.B = B
        self.LO = LO
        self.ps_length = ps_length
        self.d0 = d0
        self.srt = srt

        if echo == "Hahn":
            self._setup_echo("Hahn Echo", tau=400)
        elif echo == "Refocused":
            self._setup_echo("Refocused Echo", tau1=200, tau2=400)
        else:
            raise ValueError(
                "Only Hahn and Refocused echo's are currently supported")
        pass

    def _setup_echo(self, echo, tau=400):
        # PulseSpel_file = "/PulseSpel/phase_set"
        # run_general(self.api,
        #             [PulseSpel_file],
        #             [echo, "BrXPhase"],
        #             {"PhaseCycle": False},
        #             {"p0": self.ps_length*2, "p1": self.ps_length, "h": 20,
        #              "n": 1, "d0": self.d0, "d1": tau1, "d2": tau2, "pg": 128,
        #              "srt": self.srt},
        #             run=False
        #             )
        tau = Parameter("tau",tau,dim=4,step=0)

        
        seq = HahnEchoSequence(
            B=self.B,LO=self.LO,reptime=3e6,averages=1,shots=10, tau=tau)
        
        seq.evolution([tau])

        self.interface.launch(seq, savename="autoTUNE", start=False, tune=False)
        # PSpel = PulseSpel(sequence, self.MPFU)
        def_file, exp_file = write_pulsespel_file(seq,False,self.MPFU)
        file_path = self.interface.temp_dir + "/autoDEER_tune"
        # PSpel.save(file_path)
        save_file(file_path +'.def',def_file)
        save_file(file_path +'.exp',exp_file)

        self.api.set_field(seq.B.value)
        self.api.set_freq(seq.LO.value)


        run_general(self.api,
                    [file_path],
                    ["auto", "auto"],
                    {"PhaseCycle": True},
                    {"d0":self.d0,"pg": 64, "srt": self.srt},
                    run=False
                    )


    def tune_phase(
            self, channel: str, target: str, tol=0.1, maxiter=30) -> float:
        """Tunes the phase of a given channel to a given target using the
        standard scipy optimisation scripts. 

        Parameters
        ----------
        channel : str
            The chosen MPFU channel. Options: ['+<x>', '-<x>', '+<y>', '-<y>']
        target : str
            The target echo position, this can either be maximising (+) or
            minimising (-) either the real (R) or imaginary (I) of the echo. 
            Options: ['R+', 'R-', 'I+', 'I-']
        tol : float, optional
            The tolerance in phase parameter, by default 0.1
        maxiter : int, optional
            The maximum number of iterations in the optimisation, by default 30

        Returns
        -------
        float
            The optimal value of the phase parameter

        """

        channel_opts = ['+<x>', '-<x>', '+<y>', '-<y>']
        phase_opts = ['R+', 'R-', 'I+', 'I-']

        if channel not in channel_opts:
            raise ValueError(f'Channel must be one of: {channel_opts}')
    
        if target not in phase_opts:
            raise ValueError(f'Phase target must be one of: {phase_opts}')
        if channel == '+<x>':
            phase_channel = 'BrXPhase'
        elif channel == '-<x>':
            phase_channel = 'BrMinXPhase'
        elif channel == '+<y>':
            phase_channel = 'BrYPhase'
        elif channel == '-<y>':
            phase_channel = 'BrMinYPhase'

        if target == 'R+':
            test_fun = lambda x: -1 * np.real(x)
        elif target == 'R-':
            test_fun = lambda x: 1 * np.real(x)
        elif target == 'I+':
            test_fun = lambda x: -1 * np.imag(x)
        elif target == 'I-':
            test_fun = lambda x: 1 * np.imag(x)

        lb = 0.0
        ub = 100.0

        def objective(x, *args):
            # x = x[0]
            self.api.hidden[phase_channel].value = x  # Set phase to value
            time.sleep(self.hardware_wait)
            self.api.run_exp()
            while self.api.is_exp_running():
                time.sleep(1)
            data = self.api.acquire_scan()
            v = data.data

            val = test_fun(np.sum(v))

            print(f'Phase Setting = {x:.1f} \t Echo Amplitude = {-1*val:.2f}')

            return val

        output = minimize_scalar(
            objective, method='bounded', bounds=[lb, ub],
            options={'xatol': tol, 'maxiter': maxiter})
        result = output.x
        print(f"Optimal Phase Setting for {phase_channel} is: {result:.1f}")
        self.api.hidden[phase_channel].value = result
        return result

    def tune_power(
            self, channel: str, tol=0.1, maxiter=30,
            bounds=[0, 100]) -> float:
        """Tunes the attenuator of a given channel to a given target using the
        standard scipy optimisation scripts. 

        Parameters
        ----------
        channel : str
            The chosen MPFU channel. Options: ['+<x>', '-<x>', '+<y>', '-<y>']
        tol : float, optional
            The tolerance in attenuator parameter, by default 0.1
        maxiter : int, optional
            The maximum number of iterations in the optimisation, by default 30

        Returns
        -------
        float
            The optimal value of the attenuator parameter

        """
        channel_opts = ['+<x>', '-<x>', '+<y>', '-<y>']
        if channel not in channel_opts:
            raise ValueError(f'Channel must be one of: {channel_opts}')
        
        if channel == '+<x>':
            atten_channel = 'BrXAmp'
        elif channel == '-<x>':
            atten_channel = 'BrMinXAmp'
        elif channel == '+<y>':
            atten_channel = 'BrYAmp'
        elif channel == '-<y>':
            atten_channel = 'BrMinYAmp'

        lb = bounds[0]
        ub = bounds[1]

        def objective(x, *args):
            self.api.hidden[atten_channel].value = x  # Set phase to value
            time.sleep(self.hardware_wait)
            self.api.run_exp()
            while self.api.is_exp_running():
                time.sleep(1)
            data = self.api.acquire_scan()
            v = data.data

            val = -1 * np.sum(np.abs(v))

            print(f'Power Setting = {x:.1f} \t Echo Amplitude = {-1*val:.2f}')

            return val

        output = minimize_scalar(
            objective, method='bounded', bounds=[lb, ub],
            options={'xatol': tol, 'maxiter': maxiter})
        result = output.x
        print(f"Optimal Power Setting for {atten_channel} is: {result:.1f}")
        self.api.hidden[atten_channel].value = result
        return result

    def tune(self, channels: dict, tol: float = 0.1,
             bounds=[0, 100]) -> None:
        """Tunes both the power and attenuation for a collection of channels.

        Parameters
        ----------
        channels : dict
            A dictionary of MPFU channels to be tunned and the associated phase
            target.\\
            Channel options = ['+<x>', '-<x>', '+<y>', '-<y>']\\
            Phase target options = ['R+', 'R-', 'I+', 'I-']\\
            E.g. {'+<x>': 'R+','-<x>': 'R-'}
        tol : float, optional
            The tolerance for all optimisations, by default 0.1
        """
        for channel in channels:
            if channel == '+<x>':
                phase_cycle = 'BrXPhase'
            elif channel == '-<x>':
                phase_cycle = 'BrMinXPhase'
            elif channel == '+<y>':
                phase_cycle = 'BrYPhase'
            elif channel == '-<y>':
                phase_cycle = 'BrMinYPhase'
            self.api.set_PulseSpel_phase_cycling(phase_cycle)
            print(f"Tuning channel: {channel}")
            self.tune_power(channel, tol=tol, bounds=bounds)
            self.tune_phase(channel, channels[channel], tol=tol)

    def calc_d0(self):
        initial_d0 = 500
        PulseSpel_file = "/PulseSpel/phase_set"
        run_general(
            self.api,
            [PulseSpel_file],
            ['Hahn Echo Trans', "BrXPhase"],
            {"PhaseCycle": False},
            {"p0": self.ps_length*2, "p1": self.ps_length, "h": 20, "n": 1,
             "d1": 400, "d0": initial_d0}
            )
        while self.api.is_exp_running():
            time.sleep(1)
        data = self.api.acquire_scan()
        max_pos = np.argmax(np.abs(data.data))
        max_time = data.axes[max_pos]
        d0 = initial_d0 + max_time
        return d0
