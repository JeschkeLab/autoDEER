import matlab.engine
from scipy.io import loadmat, savemat
from autoDeer.hardware.openepr import Sequence, Pulse, RectPulse, \
    ChirpPulse, HSPulse, Detection
import numpy as np


class ETH_awg:

    def __init__(self, awg_freq=1.5, dig_rate=2) -> None:
        """Connecting to a Andrin Doll style spectrometer, commonly in use at
        ETH Zürich.

        System Requirments
        -------------------
        - Matlab 2022b or later
        - Matlab engine for python
        - Python (3.10+ recommended)
        """
        self.connect()
        self.awg_freq = awg_freq
        self.dig_rate = dig_rate
        pass

    def connect(self, session=None):

        if session is None:
            matlab_sessions = matlab.engine.find_matlab()

            if len(matlab_sessions) > 1:
                print("More than one session, picking the first.")
            session = matlab_sessions[0]
        
        self.engine = matlab.engine.connect_matlab(session)
        self.workspace = self.engine.workspace

    def acquire_dataset(self) -> dict:
        self.engine.dig_interface('savenow')
        path = ""
        loadmat(path, squeeze_me=True)

    def launch(self, experiment):
        self.workspace['exp'] = experiment
        self.engine.launch('exp', nargout=0)

    def _build_exp_struct(self, sequence) -> dict:

        struc = {}

        struc["LO"] = sequence.LO.value - self.awg_freq
        struc["avgs"] = sequence.averages.value
        struc["reptime"] = sequence.reptime.value * 1e3
        struc["shots"] = sequence.shots.value

        # Build pulse/detection events
        struc["events"] = list(map(self._build_pulse, sequence.pulses))

        # Build parvars
        struc["parvars"] = []
        struc["parvars"][0] = self._build_phase_cycle()

    def _build_pulse(self, pulse) -> dict:

        event = {}
        event["t"] = pulse.t.value

        if type(pulse) is Detection:
            event["det_pos"] = pulse.t.value
            event["det_len"] = pulse.tp.value * self.dig_rate
            pass

        # Assuming we now have an actual pulse not detection event
        event["pulsedef"] = {}
        event["pulsedef"]["scale"] = pulse.scale.value
        event["pulsedef"]["tp"] = pulse.tp.value

        if type(pulse) is RectPulse:
            event["pulsedef"]["type"] = 'chirp'
            event["pulsedef"]["nu_init"] = pulse.freq.value + self.awg_freq
        
        elif type(pulse) is ChirpPulse:
            event["pulsedef"]["type"] = 'chirp'
            
            if hasattr(pulse, "init_freq"):
                event["pulsedef"]["nu_init"] = pulse.init_freq.value +\
                     self.awg_freq
            else:
                nu_init = pulse.final_freq.value - pulse.BW.value
                event["pulsedef"]["nu_init"] = nu_init + self.awg_freq
            
            if hasattr(pulse, "final_freq"):
                event["pulsedef"]["nu_final"] = pulse.final_freq.value +\
                     self.awg_freq
            else:
                nu_final = pulse.init_freq.value + pulse.BW.value
                event["pulsedef"]["nu_final"] = nu_final + self.awg_freq
            
        elif type(pulse) is HSPulse:
            event["pulsedef"]["type"] = 'HS'
            raise RuntimeError("Not yet implemented")
        elif type(pulse) is Pulse:
            event["pulsedef"]["type"] = 'FMAM'
            raise RuntimeError("Not yet implemented")

        return event

    def _build_phase_cycle(self, sequence) -> dict:

        parvar = {}
        parvar["reduce"] = 1
        parvar["ax_id"] = 0

        pulse_str = lambda x: f"events{{{x}}}.pulsedef.phase"
        parvar["variables"] = list(map(pulse_str, sequence.pcyc_vars))

        det_id = 0
        det_str = "events{{{}}}.det_sign".format(det_id)
        parvar["variables"].append(det_str)

        parvar["vec"] = np.vstack([sequence.pcyc_cycles, sequence.pcyc_dets]).T

        return parvar

    def _build_parvar(self, id, sequence) -> dict:

        prog_table = sequence.progTable
        prog_table_n = len(prog_table[0])
        parvar = {}

        parvar["variables"] = []
        parvar["vec"] = []

        for i in range(0, prog_table_n):

            if prog_table[0][i] == id:
                pulse_num = prog_table[1][i]
                var = prog_table[2][i]
                vec = prog_table[3][i]
                if pulse_num is not None:
                    if var == "t":
                        pulse_str = f"events{{{pulse_num}}}.t"
                    else:
                        pulse_str = f"events{{{pulse_num}}}.pulsedef.{var}"
                    
                else:
                    pulse_str = var
                
                parvar["variables"].append(pulse_str)
                parvar["vec"].append(vec)

        parvar["vec"] = np.stack(parvar["vec"]).T
        return parvar
