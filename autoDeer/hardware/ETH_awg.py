import matlab.engine
from scipy.io import loadmat
from autoDeer.hardware.openepr import Sequence, Pulse, RectPulse, \
    ChirpPulse, HSPulse, Detection
import numpy as np
import os
import re


class ETH_awg_interface:

    def __init__(self, awg_freq=1.5, dig_rate=2) -> None:
        """An interface for connecting to a Andrin Doll style spectrometer,
        commonly in use at ETH Zürich.

        System Requirments
        -------------------
        - Matlab 2022b or later
        - Matlab engine for python
        - Python (3.10+ recommended)

        Parameters
        -----------
        awg_freq : float
            The normal operating AWG frequency. 
            Sequence.LO = AWG.LO + AWG.awg_freq 

        dig_rate : float
            The speed of the digitser in GSa/s
        """
        
        self.connect()
        self.awg_freq = awg_freq
        self.dig_rate = dig_rate
        pass

    def connect(self, session=None):
        """Connect to a running matlab session. If more than one session has 
        been started this will choose the first one. It is recomended that only
        one session is open at one time, or that the engine is started with a
        known name.

        Parameters
        ----------
        session : str, optional
            The string denoting a specific session to connect to
            , by default None
        """

        if session is None:
            matlab_sessions = matlab.engine.find_matlab()

            if len(matlab_sessions) > 1:
                print("More than one session, picking the first.")
            session = matlab_sessions[0]
        
        self.engine = matlab.engine.connect_matlab(session)
        self.workspace = self.engine.workspace

    def acquire_dataset(self) -> dict:
        """ Acquire and return the current or most recent dataset.

        Returns
        -------
        dict
            The dataset
        """
        cur_exp = self.workspace['currexp']
        folder_path = cur_exp['savepath']
        filename = cur_exp['savename']
        files = os.listdir(folder_path)

        def extract_date_time(str):
            output = re.findall(r"(\d{8})_(\d{4})", str)
            if output != []:
                date = int(output[0][0])
                time = int(output[0][1])
                return date*10000 + time
            else:
                return 0
        
        newest = max([extract_date_time(file) for file in files])
        date = newest // 10000
        time = newest - date * 10000
        path = folder_path + "\\" + f"{date:08d}_{time:04d}_{filename}.mat"
        
        self.engine.dig_interface('savenow')
        return loadmat(path, simplify_cells=True, squeeze_me=True)

    def launch(self, sequence: Sequence, savename: str, IFgain: int = 0):
        """Launch a sequence on the spectrometer.

        Parameters
        ----------
        sequence : Sequence
            The pulse sequence to launched.
        savename : str
            The save name for the file.
        IFgain : int
            The IF gain, either [0,1,2], default 0.
        """
        struct = self._build_exp_struct(sequence)
        struct['savename'] = savename
        struct['IFgain'] = IFgain
        self.cur_exp = struct
        self.workspace['exp'] = self.cur_exp
        self.engine.eval('launch(exp)', nargout=0)

    def _build_exp_struct(self, sequence) -> dict:

        struc = {}

        struc["LO"] = round(float(sequence.LO.value - self.awg_freq), 3)
        struc["avgs"] = float(sequence.averages.value)
        struc["reptime"] = round(float(sequence.reptime.value * 1e3), 0)
        struc["shots"] = float(sequence.shots.value)
        struc['B'] = round(float(sequence.B.value), 0)
        struc['name'] = sequence.name
        # Build pulse/detection events
        struc["events"] = list(map(self._build_pulse, sequence.pulses))

        unique_parvars = np.unique(sequence.progTable[0])

        # Build parvars
        struc["parvars"] = []
        if hasattr(sequence, 'pcyc_vars'):
            struc["parvars"].append(self._build_phase_cycle(sequence))
        for i in unique_parvars:
            struc["parvars"].append(self._build_parvar(i, sequence))
        
        return struc

    def _build_pulse(self, pulse) -> dict:

        event = {}
        event["t"] = float(pulse.t.value)

        if type(pulse) is Detection:
            event["det_len"] = float(pulse.tp.value * self.dig_rate)
            event["name"] = "det"
            return event

        # Assuming we now have an actual pulse not detection event
        event["pulsedef"] = {}
        event["pulsedef"]["scale"] = float(pulse.scale.value)
        event["pulsedef"]["tp"] = float(pulse.tp.value)

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
        parvar["ax_id"] = 1
        parvar["name"] = "phase_cycle"

        pulse_str = lambda x: f"events{{{x+1}}}.pulsedef.phase"
        parvar["variables"] = list(map(pulse_str, sequence.pcyc_vars))

        # Find detection pulse
        for i, pulse in enumerate(sequence.pulses):
            if type(pulse) == Detection:
                det_id = i

        det_str = "events{{{}}}.det_sign".format(det_id+1)
        parvar["variables"].append(det_str)

        parvar["vec"] = np.vstack([sequence.pcyc_cycles, sequence.pcyc_dets]).T

        return parvar

    def _build_parvar(self, id, sequence) -> dict:

        prog_table = sequence.progTable
        prog_table_n = len(prog_table[0])
        parvar = {}
        parvar["name"] = f"parvar{id+1}"

        parvar["variables"] = []
        parvar["vec"] = []

        for i in range(0, prog_table_n):

            if prog_table[0][i] == id:
                pulse_num = prog_table[1][i]
                var = prog_table[2][i]
                vec = prog_table[3][i]
                if pulse_num is not None:
                    if var in ["freq", "init_freq"]:
                        vec += self.awg_freq
                        var = "nu_init"
                    if var == "final_freq":
                        vec += self.awg_freq
                        var = "nu_final"
            
                    if var == "t":
                        pulse_str = f"events{{{pulse_num+1}}}.t"
                    else:
                        pulse_str = f"events{{{pulse_num+1}}}.pulsedef.{var}"

                else:
                    pulse_str = var
                
                parvar["variables"].append(pulse_str)
                parvar["vec"].append(vec)

        parvar["vec"] = np.stack(parvar["vec"]).T
        return parvar
