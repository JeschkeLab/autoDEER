import matlab.engine
from autodeer.dataset import Dataset
from autodeer.classes import  Interface, Parameter
from autodeer.pulses import Pulse, RectPulse, ChirpPulse, HSPulse, Delay, Detection
from autodeer.sequences import Sequence, HahnEchoSequence
import numpy as np
import os
import re
import time
from deerlab import correctphase
import numpy as np
import scipy.signal as sig
from scipy.io import loadmat
from autodeer.utils import transpose_list_of_dicts


class ETH_awg_interface(Interface):
    """
    Represents the interface for connecting to Andrin Doll style spectrometers.
    """
    def __init__(self, awg_freq=1.5, dig_rate=2, test_mode=False) -> None:
        """An interface for connecting to a Andrin Doll style spectrometer,
        commonly in use at ETH Zürich.

        System Requirements
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
        if not test_mode:
            self.connect()
            
        self.awg_freq = awg_freq
        self.dig_rate = dig_rate
        self.pulses = {}
        self.cur_exp = None
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
            elif len(matlab_sessions) == 0:
                raise RuntimeError(
                    "A matlab python session must be started. \n"+
                    "Please type into matlab session: "
                    "matlab.engine.shareEngine"
                    )
        
        self.engine = matlab.engine.connect_matlab(session)
        self.workspace = self.engine.workspace

    def acquire_dataset(self, options={}, verbosity=0) -> Dataset:
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
                start_time = int(output[0][1])
                return date*10000 + start_time
            else:
                return 0
        
        newest = max([extract_date_time(file) for file in files])
        date = newest // 10000
        start_time = newest - date * 10000
        path = folder_path + "\\" \
            + f"{date:08d}_{start_time:04d}_{filename}.mat"
        
        self.engine.dig_interface('savenow')
        
        for i in range(0, 10):
            try:
                Matfile = loadmat(path, simplify_cells=True, squeeze_me=True)
                data = uwb_load(Matfile, options=options, verbosity=verbosity)
                data.LO = Parameter("LO", data.params['LO']+self.awg_freq, unit="GHz", description="Total local oscilator frequency")
                data.sequence = self.cur_exp
            except OSError:
                time.sleep(10)
            else:
                return data
        
    def launch(self, sequence , savename: str, IFgain: int = 0):
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
        self.cur_exp = sequence
        self.workspace['exp'] = struct
        self.engine.eval('launch(exp)', nargout=0)

    def isrunning(self) -> bool:
        state = bool(self.engine.dig_interface('savenow'))
        return state
    
    def tune_rectpulse(self,*,tp, LO, B, reptime):
        """Generates a rectangular pi and pi/2 pulse of the given length at 
        the given field position. This value is stored in the pulse cache. 

        Parameters
        ----------
        tp : float
            Pulse length in ns
        LO : float
            Central frequency of this pulse in GHz
        B : float
            Magnetic B0 field position in Gauss
        reptime: float
            Shot repetion time in us.

        Returns
        -------
        p90: RectPulse
            A tuned rectangular pi/2 pulse of length tp
        p180: RectPulse
            A tuned rectangular pi pulse of length tp
        """

        amp_tune =HahnEchoSequence(
            B=B, LO=LO, reptime=reptime, averages=1, shots=400
        )

        scale = Parameter("scale",0,dim=45,step=0.02)
        amp_tune.pulses[0].tp.value = tp
        amp_tune.pulses[0].scale = scale
        amp_tune.pulses[1].tp.value = tp * 2
        amp_tune.pulses[1].scale = scale

        amp_tune.evolution([scale])
        
        # amp_tune.addPulsesProg(
        #     pulses=[0,1],
        #     variables=['scale','scale'],
        #     axis_id=0,
        #     axis=np.arange(0,0.9,0.02),
        # )

        self.launch(amp_tune, "autoDEER_amptune", IFgain=1)

        while self.isrunning():
            time.sleep(10)
        dataset = self.acquire_dataset()
        scale = np.around(dataset.axes[0][dataset.data.argmax()],2)
        if scale > 0.9:
            raise RuntimeError("Not enough power avaliable.")
        
        self.pulses[f"p90_{tp}"] = amp_tune.pulses[0].copy(
            scale=scale, LO=amp_tune.LO)
        self.pulses[f"p180_{tp*2}"] = amp_tune.pulses[1].copy(
            scale=scale, LO=amp_tune.LO)

        return self.pulses[f"p90_{tp}"], self.pulses[f"p180_{tp*2}"]

    
    def tune_pulse(self, pulse, mode, LO, B , reptime):
        """Tunes a single pulse a range of methods.

        Parameters
        ----------
        pulse : Pulse
            The Pulse object in need of tuning.
        mode : str
            The method to be used.
        LO : float
            The local oscilator frequency in GHz
        B : float
            Magnetic B0 field position in Gauss
        reptime : us
            Shot repetion time in us.

        Returns
        -------
        Tunned Pulse: Pulse
            The returned pulse object that is now tunned.
        """
        # Check pulse is a pulse
        if type(pulse) == Delay:
            pass
        if type(pulse) == Detection:
            pass
        
        # Get absolute central frequency
        if hasattr(pulse,"freq"):
            c_frq = pulse.freq.value + LO
        elif hasattr(pulse, "init_freq") & hasattr(pulse, "BW"):
            c_frq = pulse.init_freq.value + 0.5*pulse.BW.value + LO
        elif hasattr(pulse, "final_freq") & hasattr(pulse, "BW"):
            c_frq = pulse.final_freq.value - 0.5*pulse.BW.value + LO
        elif hasattr(pulse, "init_freq") & hasattr(pulse, "final_freq"):
            c_frq = 0.5*(pulse.final_freq.value + pulse.final_freq.value) + LO

        # Find rect pulses

        all_rect_pulses = list(self.pulses.keys())
        pulse_matches = []
        for pulse_name in all_rect_pulses:
            # if not re.match(r"^p180_",pulse_name):
            #     continue
            pulse_frq = self.pulses[pulse_name].LO.value + self.pulses[pulse_name].freq.value
            if not np.abs(pulse_frq - c_frq) < 0.01: #Within 10MHz
                continue
            pulse_matches.append(pulse_name)
        
        # Find best pi pulse
        pi_length_best = 1e6
        for pulse_name in pulse_matches:
            if re.match(r"^p180_",pulse_name):
                ps_length = int(re.search(r"p180_(\d+)",pulse_name).groups()[0])
                if ps_length < pi_length_best:
                    pi_length_best = ps_length
        
        if pi_length_best == 1e6:
            _, pi_pulse = self.tune_rectpulse(tp=12, B=B, LO=LO, reptime=reptime)
        else:
            pi_pulse = self.pulses[f"p180_{pi_length_best}"]

        pi2_length_best = 1e6
        for pulse_name in pulse_matches:
            if re.match(r"^p90_",pulse_name):
                ps_length = int(re.search(r"p90_(\d+)",pulse_name).groups()[0])
                if ps_length < pi2_length_best:
                    pi2_length_best = ps_length
        
        if pi2_length_best == 1e6:
            pi2_pulse, _ = self.tune_rectpulse(tp=12, B=B, LO=LO,reptime=reptime)
        else:
            pi2_pulse = self.pulses[f"p90_{pi2_length_best}"]


        if mode == "amp_hahn":
            amp_tune =HahnEchoSequence(
                B=B, LO=LO, 
                reptime=reptime, averages=1, shots=400,
                pi2_pulse = pulse, pi_pulse=pi_pulse
            )

            scale = Parameter('scale',0,unit=None,step=0.02, dim=45, description='The amplitude of the pulse 0-1')
            amp_tune.pulses[0].scale = scale

            amp_tune.evolution([scale])

            self.launch(amp_tune, "autoDEER_amptune", IFgain=1)

            while self.isrunning():
                time.sleep(10)
            dataset = self.acquire_dataset()
            new_amp = np.around(dataset.axes[0][dataset.data.argmax()],2)
            pulse.scale = Parameter('scale',new_amp,unit=None,description='The amplitude of the pulse 0-1')
            return pulse

        elif mode == "amp_nut":

            nut_tune = Sequence(
                name="nut_tune", B=B, LO=LO, reptime=reptime,
                averages=1,shots=400
            )
            nut_tune.addPulse(pulse.copy(
                t=0, pcyc={"phases":[0],"dets":[1]}, scale=0))
            nut_tune.addPulse(
                pi2_pulse.copy(t=2e3,
                               pcyc={"phases":[0, np.pi],"dets":[1, -1]},
                               freq=c_frq-LO))
            nut_tune.addPulse(
                pi_pulse.copy(t=2.5e3, pcyc={"phases":[0],"dets":[1]},
                              freq=c_frq-LO))
            nut_tune.addPulse(Detection(t=3e3, tp=512, freq=c_frq-LO))

            scale = Parameter('scale',0,unit=None,step=0.02, dim=45, description='The amplitude of the pulse 0-1')
            nut_tune.pulses[0].scale = scale
            nut_tune.evolution([scale])


            # nut_tune.addPulsesProg(
            #     pulses=[0],
            #     variables=["scale"],
            #     axis_id = 0,
            #     axis= np.arange(0,0.9,0.02)
            # )
            self.launch(nut_tune, "autoDEER_amptune", IFgain=1)

            while self.isrunning():
                time.sleep(10)
            dataset = self.acquire_dataset()
            data = correctphase(dataset.data)
            if data[0] < 0:
                data *= -1

            if np.isclose(pulse.flipangle.value, np.pi):
                new_amp = np.around(dataset.axes[0][data.argmin()],2)
            elif np.isclose(pulse.flipangle.value, np.pi/2):
                sign_changes = np.diff(np.sign(np.real(data)))
                new_amp = np.around(dataset.axes[0][np.nonzero(sign_changes)[0][0]],2)
            else:
                raise RuntimeError("Target pulse can only have a flip angle of either: ",
                                "pi or pi/2.")
            pulse.scale = Parameter('scale',new_amp,unit=None,description='The amplitude of the pulse 0-1')
        
            return pulse
    
    def tune(self,*, sequence=None, mode="amp_hahn", LO=None, gyro=None):

        if mode == "rect_tune":
            if LO is None:
                raise ValueError("LO must be given for rect_tune")
            if gyro is None:
                raise ValueError("gyro must be give")
            elif gyro >1:
                raise ValueError("Gyromagnetic ratio must be give in GHz/G")
            
            amp_tune =HahnEchoSequence(
                B=LO/gyro, LO=LO, reptime=2e3, averages=1, shots=400
            )
            tp = 12
            amp_tune.pulses[0].tp.value = tp
            amp_tune.pulses[0].scale.value = 0
            amp_tune.pulses[1].tp.value = tp*2
            amp_tune.pulses[1].scale.value = 0
            
            amp_tune.addPulsesProg(
                pulses=[0,1],
                variables=['scale','scale'],
                axis_id=0,
                axis=np.arange(0,0.9,0.02),
            )

            self.launch(amp_tune, "autoDEER_amptune", IFgain=1)

            while self.isrunning():
                time.sleep(10)
            dataset = self.acquire_dataset()
            scale = np.around(dataset.axes[0][dataset.data.argmax()],2)
            if scale > 0.9:
                raise RuntimeError("Not enough power avaliable.")
            
            self.pulses[f"p90_{tp}"] = amp_tune.pulses[0].copy(
                scale=scale, LO=amp_tune.LO)
            self.pulses[f"p180_{tp*2}"] = amp_tune.pulses[1].copy(
                scale=scale, LO=amp_tune.LO)
        
        elif mode == "amp_hahn":
            for pulse in sequence.pulses:
                if type(pulse) == Delay:
                    continue
                if type(pulse) == Detection:
                    continue

                all_pulses = list(self.pulses.keys())
                pulse_matches = []
                for pulse_name in all_pulses:
                    if not re.match(r"^p180_",pulse_name):
                        continue
                    if not np.abs((self.pulses[pulse_name].LO.value + self.pulses[pulse_name].freq.value) - (sequence.LO.value + pulse.freq.value)) < 0.01:
                        continue
                    pulse_matches.append(pulse_name)
                    
                    ps_length_best =1e6
                for pulse_name in pulse_matches:
                    ps_length = int(re.search(r"p180_(\d+)",pulse_name).groups()[0])
                    if ps_length < ps_length_best:
                        ps_length_best = ps_length
                
                pi_pulse = self.pulses[f"p180_{ps_length_best}"]
                

                amp_tune =HahnEchoSequence(
                    B=sequence.B.value, LO=sequence.LO.value, 
                    reptime=sequence.reptime.value, averages=1, shots=400,
                    pi2_pulse = pulse, pi_pulse=pi_pulse
                )
        
                amp_tune.pulses[0].scale.value = 0

                amp_tune.addPulsesProg(
                    pulses=[0],
                    variables=['scale'],
                    axis_id=0,
                    axis=np.arange(0,0.9,0.02),
                )

                self.launch(amp_tune, "autoDEER_amptune", IFgain=1)

                while self.isrunning():
                    time.sleep(10)
                dataset = self.acquire_dataset()
                scale = np.around(dataset.axes[0][dataset.data.argmax()],2)
                pulse.scale.value = scale

            return sequence
                    

    def _build_exp_struct(self, sequence) -> dict:

        struc = {}

        struc["avgs"] = float(sequence.averages.value)
        struc["reptime"] = round(float(sequence.reptime.value * 1e3), 0)
        struc["shots"] = float(sequence.shots.value)
        struc['B'] = round(float(sequence.B.value), 0)
        struc['name'] = sequence.name
        # Build pulse/detection events
        struc["events"] = list(map(self._build_pulse, sequence.pulses))

        unique_parvars = np.unique(sequence.progTable["axID"])

        # Build parvars
        struc["parvars"] = []
        if hasattr(sequence, 'pcyc_vars'):
            pcyc_parvar = self._build_phase_cycle(sequence)
            if pcyc_parvar["vec"].shape[0] != 1:
                struc["parvars"].append(pcyc_parvar)
        for i in unique_parvars:
            struc["parvars"].append(self._build_parvar(i, sequence))
        
        struc["LO"] = round(float(sequence.LO.value - self.awg_freq), 3)

        return struc

    def _build_pulse(self, pulse) -> dict:

        event = {}
        event["t"] = float(pulse.t.value)

        if type(pulse) is Detection:
            event["det_len"] = float(pulse.tp.value * self.dig_rate)
            event["det_frq"] = float(pulse.freq.value) + self.awg_freq
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

            event["pulsedef"]["HSorder"] = float(pulse.order1.value)
            event["pulsedef"]["HSorder2"] = float(pulse.order2.value)
            event["pulsedef"]["HSbeta"] = float(pulse.beta.value)

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

        parvar["vec"] = np.vstack(
            [sequence.pcyc_cycles.T, sequence.pcyc_dets]).T

        return parvar

    def _build_parvar(self, id, sequence) -> dict:
        """This interface takes a dictionary called a parvar for all 
        progressive elements. It is this object that controls how the
        sequence changes with time.


        .. note::
            This interface interprets any change in `LO` as being a 
            change in the IF frequency of all pulses and detections.
             I.e. the physcial LO **does not** change. 

        Parameters
        ----------
        id : _type_
            _description_
        sequence : _type_
            _description_

        Returns
        -------
        dict
            _description_
        """

        def get_vec(seq,EventID,variable,uuid):

            if EventID is None:
                attr = getattr(seq, variable)
            else:
                attr = getattr(seq.pulses[EventID], variable)

            # Find local position
            axis_T = transpose_list_of_dicts(attr.axis)
            i = axis_T['uuid'].index(uuid)
            vec = attr.value + attr.axis[i]['axis']
            return vec.astype(float)

        prog_table = sequence.progTable
        prog_table_n = len(prog_table["axID"])
        parvar = {}
        parvar["name"] = f"parvar{id+1}"

        parvar["variables"] = []
        parvar["vec"] = []

        for i in range(0, prog_table_n):

            if prog_table["axID"][i] == id:
                pulse_num = prog_table["EventID"][i]
                var = prog_table["Variable"][i]
                uuid = prog_table["uuid"][i]
                # vec = prog_table["axis"][i].astype(float)
                vec = get_vec(sequence,pulse_num,var,uuid)
                if pulse_num is not None:
                    if var in ["freq", "init_freq"]:
                        vec += self.awg_freq
                        var = "nu_init"
                    elif var == "final_freq":
                        vec += self.awg_freq
                        var = "nu_final"

                    if var == "t":
                        pulse_str = f"events{{{pulse_num+1}}}.t"
                    else:
                        pulse_str = f"events{{{pulse_num+1}}}.pulsedef.{var}"
                    
                    parvar["variables"].append(pulse_str)
                    parvar["vec"].append(vec)

                elif var == "LO":
                    # Instead of stepping the LO we will step the frequency
                    # all the pulses and detection events.
                    pulse_strings = []
                    # vec = prog_table["axis"][i].astype(float)
                    centre_freq = (vec[-1] + vec[0])/2
                    LO = centre_freq - self.awg_freq
                    sequence.LO.value = centre_freq
                    vec = vec - LO
                    vecs = []
                    pulse_str = lambda i: f"events{{{i+1}}}.pulsedef.nu_init"
                    det_str = lambda i: f"events{{{i+1}}}.det_frq"
                    for i,pulses in enumerate(sequence.pulses):
                        if type(pulses) is Delay:
                            continue
                        elif type(pulses) is Detection:
                            pulse_strings.append(det_str(i))
                            vecs.append(vec)
                        else:
                            pulse_strings.append(pulse_str(i))
                            vecs.append(vec)

                    parvar["variables"].extend(pulse_strings)
                    parvar["vec"].extend(vecs)                
                
                else:
                    pulse_str = var
                    parvar["variables"].append(pulse_str)
                    parvar["vec"].append(vec)
                


        parvar["vec"] = np.stack(parvar["vec"]).T
        return parvar

    def terminate(self):
        """ Stops the current experiment
        """
        self.engine.dig_interface('terminate', nargout=0)
        pass


def uwb_load(matfile: np.ndarray, options: dict = dict(), verbosity=0):
    """uwb_load This function is based upon the uwb_eval matlab function 
    developed by Andrin Doll. This is mostly a Python version of this software,
     for loading data generated by Doll style spectrometers.

    Parameters
    ----------
    matfile : np.ndarray
        _description_
    options : dict, optional
        _description_, by default None
    """

    # Extract Data
    estr = matfile[matfile['expname']]
    conf = matfile['conf'] 
    
    def extract_data(matfile):
        if "dta" in matfile.keys():
            nAvgs = matfile["nAvgs"]
            dta = [matfile["dta"]]
        elif "dta_001" in matfile.keys():
            dta = []
            for ii in range(1, estr["avgs"]+1):
                actname = 'dta_%03u' % ii 
                if actname in matfile.keys():
                    dta.append(matfile[actname])
                    # Only keep it if the average is complete, unless it is 
                    # the first
                    if sum(dta[ii-1][:, -1]) == 0 and ii > 1:
                        dta = dta[:-1]
                    elif sum(dta[ii-1][:, -1]) == 0 and ii == 1:
                        nAvgs = 0
                    else:
                        nAvgs = ii
        return [dta, nAvgs]
    
    dta, nAvgs = extract_data(matfile)

    # Eliminate Phase cycles

    if "postsigns" not in estr.keys():
        print("TODO: check uwb_eval")
        raise RuntimeError("This is not implemented yet")

    if np.isscalar(estr["postsigns"]["signs"]):
        estr["postsigns"]["signs"] = [estr["postsigns"]["signs"]]

    if type(estr["parvars"]) is dict:
        estr["parvars"] = [estr["parvars"]]

    cycled = np.array(list(map(np.size, estr["postsigns"]["signs"]))) > 1

    #  decide on wheteher the phase cycle should be eliminated or not
    if any(cycled == 0):
        elim_pcyc = 1  # if there is any non-phasecycling parvar
    else:
        elim_pcyc = 0  # if all parvars cycle phases, do not reduce them 
    
    if "elim_pcyc" in options.keys():
        elim_pcyc = options["elim_pcyc"]
    
    # Get the cycles out
    if elim_pcyc:
        for ii in range(0, len(cycled)):
            if cycled[ii]:
                if ii > 1:
                    n_skip = np.prod(estr["postsigns"]["dims"][0:ii-1])
                else:
                    n_skip = 1
                plus_idx = np.where(estr["postsigns"]["signs"][ii] == 1)
                minus_idx = np.where(estr["postsigns"]["signs"][ii] == -1)
                plus_mask = np.arange(0, n_skip) + (plus_idx - 1) * n_skip
                minus_mask = np.arange(0, n_skip) + (minus_idx - 1) * n_skip
                n_rep = np.size(dta[1], 2) / \
                    (n_skip * estr["postsigns"]["dims"][ii])

                for kk in range(0, len(dta)):
                    #  re-allocate
                    tmp = dta[kk]
                    dta[kk] = np.zeros(np.size(tmp, 0), n_rep*n_skip)
                    # subtract out
                    for jj in range(0, n_rep):
                        curr_offset = (jj) * n_skip
                        full_offset = (jj) * n_skip * \
                            estr["postsigns"]["dims"][ii]
                        dta[kk][:, np.arange(0, n_skip)+curr_offset] = \
                            tmp[:, plus_mask+full_offset] - \
                            tmp[:, minus_mask+full_offset]

    #  Find all the axes
    dta_x = []
    ii_dtax = 0
    relevant_parvars = []

    for ii in range(0, len(cycled)):
        estr["postsigns"]["ids"] = \
            np.array(estr["postsigns"]["ids"]).reshape(-1)
        if not (elim_pcyc and cycled[ii]):
            if type(estr["parvars"]) is list:
                dta_x.append(estr["parvars"][estr["postsigns"]["ids"][ii]-1]
                            ["axis"].astype(np.float32))
            else:
                dta_x.append(estr["parvars"]["axis"].astype(np.float32))
            relevant_parvars.append(estr["postsigns"]["ids"][ii]-1)
            ii_dtax += 1
    
    exp_dim = ii_dtax
    if ii_dtax == 0:
        raise RuntimeError("Your dataset does not have any swept dimensions." 
                           "Uwb_eval does not work for experiments without"
                           "any parvars")
    elif ii_dtax > 2:
        raise RuntimeError("Uwb_eval cannot handle more than two dimensions")
    
    det_frq = estr["events"][estr["det_event"]-1]["det_frq"]
    det_frq_dim = 0
    fsmp = conf["std"]["dig_rate"]

    # Check for any frequency changes as well as any fixed downconversion 
    # frequencies

    if "det_frq" in options.keys():
        det_frq = options["det_frq"]
    else:

        det_frq_dim = 0
        for ii in range(0, len(relevant_parvars)):
            act_par = estr["parvars"][relevant_parvars[ii]]
            frq_change = np.zeros((len(act_par["variables"]), 1))
            for jj in range(0, len(act_par["variables"])):
                if not any('nu_' in word for word
                           in estr["parvars"][0]["variables"]):
                    frq_change[jj] = 1

        if any(frq_change):  # was there a frequency change
            # is the frequency change relevant
            if "det_frq_id" in estr["events"][estr["det_event"]-1]:
                frq_pulse = estr["events"][estr["det_event"]-1]["det_frq_id"]
                nu_init_change = 0
                nu_final_change = 0
                for jj in range(0, np.size(act_par["variables"])):
                    if ('events{' + str(frq_pulse+1) + "}.pulsedef.nu_i") in \
                            act_par["variables"][jj]:
                        nu_init_change = jj
                    if ('events{' + str(frq_pulse+1) + "}.pulsedef.nu_f") in \
                            act_par["variables"][jj]:
                        nu_final_change = jj

                if any([nu_init_change, nu_final_change]):
                    # There is a frequency change on the frequency encoding 
                    # pulse

                    # dimension that will determine the detection frequency
                    det_frq_dim = ii  

                    if "nu_final" != estr["events"][frq_pulse]['pulsedef']:
                        # rectangular pulse
                        
                        if nu_init_change == 0:
                            print(
                                "uwb_eval has no idea how to guess your "
                                "detection frequency. You were setting a "
                                "rectangular pulse in event" + str(frq_pulse) 
                                + ", but are now increasing its end frequency"
                                ". You may obtain unexpected results.")

                        # get the frequencies either from 
                        # the vectorial definition
                        if "vec" in act_par.keys():
                            det_frq = act_par["vec"][:, nu_init_change]
                        else:
                            # or from the parametric definition
                            nu_init = estr["events"][frq_pulse]["pulsedef"][
                                "nu_init"]
                            if np.isnan(act_par["strt"][nu_init_change]):
                                det_frq = np.arange(0, act_par["dim"]) * \
                                    act_par["inc"][nu_init_change] + nu_init
                            else:
                                det_frq = np.arange(0, act_par["dim"]) * \
                                    act_par["inc"][nu_init_change] + \
                                    act_par["strt"][nu_init_change]
                    
                    else:
                        #  chirp pulse, both nu_init and nu_final need to be 
                        # considered

                        nu_init = estr["events"][frq_pulse]["pulsedef"][
                            "nu_init"]
                        nu_final = estr["events"][frq_pulse]["pulsedef"][
                            "nu_final"]

                        # get the frequencies either from the 
                        # vectorial definition
                        if "vec" in act_par.keys():
                            if nu_init_change != 0:
                                nu_init = act_par["vec"][:, nu_init_change]
                            if nu_final_change != 0:
                                nu_final = act_par["vec"][:, nu_final_change]
                        else:
                            # or from the parametric definition
                            if nu_init_change != 0:
                                if np.isnan(act_par["strt"][nu_init_change]):
                                    nu_init = np.arange(0, act_par["dim"]) * \
                                        act_par["inc"][nu_init_change] + \
                                        nu_init
                                else:
                                    nu_init = np.arange(0, act_par["dim"]) * \
                                        act_par["inc"][nu_init_change] + \
                                        act_par["strt"][nu_init_change]
                            if nu_final_change != 0:
                                if np.isnan(act_par["strt"][nu_init_change]):
                                    nu_final = np.arange(0, act_par["dim"]) * \
                                        act_par["inc"][nu_final_change] + \
                                        nu_final
                                else:
                                    nu_final = np.arange(0, act_par["dim"]) * \
                                        act_par["inc"][nu_final_change] + \
                                        act_par["strt"][nu_final]
                        
                        det_frq = (nu_init + nu_final) / 2
            else:
                # we can only land here, if there was no det_frq_id given, but
                # det_frq was explicitly provided in the experiment. This could
                # be intentional, but could even so be a mistake of the user.
                print("uwb_eval has no idea how to guess your detection"
                      "frequency. You were changing some pulse frequencies, "
                      "but did not provide det_frq_id for your detection event"
                      ". I will use det_frq, as you provided it in the "
                      "experiment.")

    #  ****Check digitizer level

    parvar_pts = np.zeros(np.size(estr["parvars"]))
    for ii in range(0, len(estr["parvars"])):
        if "vec" in estr["parvars"][ii]:
            parvar_pts[ii] = np.size(estr["parvars"][ii]["vec"], 0)
        else:
            parvar_pts[ii] = estr["parvars"][ii]["dim"]
    
    # the number of data points entering one echo transient (due to reduction
    # during acquisition or reduction of phasecycles just above)
    n_traces = np.prod(parvar_pts) / np.prod(list(map(len, dta_x)))

    if "dig_max" in conf.keys():
        trace_maxlev = n_traces * estr["shots"] * conf["dig_max"]
    else:
        trace_maxlev = n_traces * estr["shots"] * 2**11
    
    #  ***** Extract all the echoes
    echopos = estr["events"][estr["det_event"]-1]["det_len"]/2 - estr[
        "events"][estr["det_event"]-1]["det_pos"] * fsmp
    dist = min([echopos, estr["events"][estr["det_event"]-1]["det_len"] - 
                echopos])
    ran_echomax = np.arange(echopos - dist, echopos + dist, dtype=np.int64)

    # get the downconversion to LO
    t_ax_full = np.arange(0, len(ran_echomax)) / fsmp
    if not np.isscalar(det_frq):
        tf = np.matmul(t_ax_full[:, None], det_frq[None, :])
        LO = np.exp(-2 * np.pi * 1j * tf)
    else:
        LO = np.exp(-2 * np.pi * 1j * t_ax_full * det_frq)

    flipback = 0

    # 1D or 2D
    if exp_dim == 2:
        dta_ev = np.zeros((len(dta_x[0]), len(dta_x[1])), dtype=np.complex128)
        dta_avg = np.zeros(
            (len(ran_echomax), len(dta_x[0]), len(dta_x[1])),
            dtype=np.complex128)
        perm_order = [0, 1, 2]
        if det_frq_dim == 1:
            perm_order = [0, 2, 1]
            flipback = 1
            dta_ev = np.transpose(dta_ev)
            dta_avg = np.transpose(dta_avg, perm_order)
    elif exp_dim == 1:
        dta_ev = np.zeros((np.size(dta_x[0])), dtype=np.complex128)
        dta_avg = np.zeros((len(ran_echomax), np.size(dta_x[0])),
                           dtype=np.complex128)


    dta_scans = np.zeros((len(dta),) + dta_ev.shape, dtype=np.complex128)
    
    for ii in range(0, len(dta)):
        dta_c = dta[ii][ran_echomax, :]
        dta_c = np.conj(np.apply_along_axis(sig.hilbert, 0, dta_c))

        # reshape the 2D data
        if exp_dim == 2:
            dta_resort = np.reshape(dta_c, (len(ran_echomax), len(dta_x[0]), 
                                    len(dta_x[1])), order='F')
            dta_resort = np.transpose(dta_resort, perm_order)
        else:
            dta_resort = dta_c

        # downconvert
        dta_dc = (dta_resort.T * LO.T).T

        # refine the echo position for further evaluation, 
        # based on the first average
        if ii == 0:

            # put a symetric window to mask the expected echo position
            window = sig.windows.chebwin(np.size(dta_dc, 0), 100)
            dta_win = np.transpose(np.transpose(dta_dc) * window)
            
            # absolute part of integral, since phase is not yet determined
            absofsum = np.squeeze(np.abs(np.sum(dta_win, 0)))

            # get the strongest echo of that series
            ref_echo = np.argmax(absofsum.flatten('F'))

            # use this echo to inform about the digitizer scale

            max_amp = np.amax(dta[ii][ran_echomax, ref_echo], 0)
            dig_level = max_amp / trace_maxlev

            if "IFgain_levels" in conf["std"]:
                # Check about eventual improvemnets by changing IF levels
                possible_levels = dig_level * conf["std"]["IFgain_levels"] / \
                    conf["std"]["IFgain_levels"][estr["IFgain"]]

                possible_levels[possible_levels > 0.75] = 0
                best_lev = np.amax(possible_levels)
                best_idx = np.argmax(possible_levels)
                if (best_idx != estr["IFgain"]) & (verbosity > 0):
                    print(
                        f"You are currently using {dig_level} of the maximum"
                        f"possible level of the digitizer at an IFgain setting"
                        f"of {estr['IFgain']} \n It may be advantageous to use"
                        f"an IFgain setting of {best_idx} , where the maximum "
                        f"level will be on the order of {best_lev}.")

            # for 2D data, only a certain slice may be requested

            if "ref_echo_2D_idx" in options.keys():
                if "ref_echo_2D_dim" not in options.keys():
                    options["ref_echo_2D_dim"] = 1
                
                if flipback:
                    if options["ref_echo_2D_dim"] == 1:
                        options["ref_echo_2D_dim"] = 2
                    else:
                        options["ref_echo_2D_dim"] = 1
                
                if options["ref_echo_2D_idx"] == "end":
                    options["ref_echo_2D_idx"] = np.size(
                        absofsum, options["ref_echo_2D_dim"]-1)

                if options["ref_echo_2D_dim"] == 1:
                    ii_ref = options["ref_echo_2D_idx"] - 1
                    jj_ref = np.argmax(absofsum[ii_ref, :])
                else:
                    jj_ref = options["ref_echo_2D_idx"] - 1
                    ii_ref = np.argmax(absofsum[jj_ref, :])
                # convert the ii_ref,jj_ref to a linear index, as this is how
                # the convolution is done a few lines below

                ref_echo = np.ravel_multi_index(
                    [ii_ref, jj_ref], absofsum.shape)

            if "ref_echo" in options.keys():
                if "end" == options["ref_echo"]:
                    ref_echo = len(absofsum[:])
                else:
                    ref_echo = options["ref_echo"]
            
            # look for zerotime by crosscorrelation with 
            # echo-like window (chebwin..),
            # use conv istead of xcorr, because of the 
            # life-easy-making 'same' option
            # TODO turn this into a matched filter
            # this is where one could use a matched echo shape 

            convshape = sig.windows.chebwin(min([100, len(ran_echomax)]), 100) 

            if exp_dim == 2:
                ref_echo_unravel = np.unravel_index(ref_echo, absofsum.shape)
                e_idx = np.argmax(sig.convolve(
                    np.abs(dta_dc[:, ref_echo_unravel[0],
                    ref_echo_unravel[1]]), convshape, mode="same"))

            else:
                e_idx = np.argmax(sig.convolve(
                    np.abs(dta_dc[:, ref_echo]), convshape, mode="same"))

            # now get the final echo window, which is centered around the
            # maximum position just found

            dist = min([e_idx, np.size(dta_dc, 0)-e_idx])
            evlen = 2 * dist

            if "evlen" in options.keys():
                evlen = options["evlen"]
            if "find_echo" in options.keys():
                e_idx = np.floor(np.size(dta_dc, 0)/2)
            
            # here the final range...
            ran_echo = np.arange(e_idx-evlen/2, e_idx+evlen/2, dtype=np.int16)
            # ... and a check wether this is applicable

            if not (ran_echo[0] >= 0 and ran_echo[-1] <= np.size(dta_dc, 0)):
                raise RuntimeError(
                    f"Echo position at {e_idx} with evaluation length of "
                    f"{evlen} is not valid, since the dataset has only "
                    f"{np.size(dta_dc,0)} points.")

            # here the final time axis of the dataset
            t_ax = np.arange(-evlen/2, evlen/2) / fsmp

            # get also indices of reference echo in case of 2D data
            if absofsum.ndim == 2:
                [ii_ref, jj_ref] = np.unravel_index(
                    ref_echo, absofsum.shape, order='F')
        
        # window the echo
        if exp_dim == 2:
            dta_win = np.multiply(dta_dc[ran_echo, :, :].T, 
                                  sig.windows.chebwin(evlen, 100)).T
        else:
            dta_win = np.multiply(dta_dc[ran_echo, :].T, 
                                  sig.windows.chebwin(evlen, 100)).T

        # get all the phases and use reference echo for normalization
        dta_ang = np.angle(np.sum(dta_win, 0))

        # for frequency changes, we phase each frequency

        if det_frq_dim != 0:
            if exp_dim == 2:
                corr_phase = dta_ang[..., jj_ref]
            else:
                corr_phase = dta_ang
        else:
            corr_phase = dta_ang[ref_echo]
        
        # check if a fixed phase was provided
        if "corr_phase" in options.keys():
            corr_phase = options["corr_phase"]
        
        # check if any datapoint needs to be phased individually
        if "phase_all" in options.keys() and options["phase_all"] == 1:
            corr_phase = dta_ang
        
        bfunc = lambda x: x * np.exp(-1j * corr_phase)
        # dta_pha = np.multiply(dta_win, np.exp(-1j * corr_phase))
        
        dta_pha = np.apply_along_axis(bfunc, 1, dta_win)
        
        dta_this_scan = np.squeeze(np.sum(dta_pha, 0)) / \
            sum(sig.windows.chebwin(evlen, 100))
        dta_ev = dta_ev + dta_this_scan

        dta_scans[ii, :] = dta_this_scan  # This will not work for 2D
        
        if exp_dim == 2:
            dta_avg[0:evlen, :, :] = dta_avg[0:evlen, :, :] + \
                np.apply_along_axis(bfunc, 1, dta_win)
#                np.multiply(dta_dc[ran_echo, :, :], np.exp(-1j * corr_phase))
        else:
            dta_avg[0:evlen, :] = dta_avg[0:evlen, :] + \
                np.multiply(dta_dc[ran_echo, :], np.exp(-1j * corr_phase))

    dta_avg = dta_avg[0:evlen, ...]
    # keyboard
    # flip back 2D Data
    if flipback:
        dta_avg = np.transpose(dta_avg, perm_order)
        dta_ev = np.transpose(dta_ev)

    output = Dataset(
        axes=dta_x,
        data=dta_ev,
        params=estr)
    
    output.scans = dta_scans
    output.add_variable(Parameter(name='nAvgs', value=nAvgs))

    # output = AWGdata(t_ax, dta_avg)
    # output.nAvgs = nAvgs
    # output.dta_x = dta_x
    # output.dta_ev = dta_ev
    # output.dta_scans = dta_scans
    # output.exp = estr
    # output.det_frq = det_frq
    # output.echopos = echopos
    # output.corr_phase = corr_phase
    # output.dig_level = dig_level
    # output.raw_dta = dta
    
    return output


# class AWGdata:

#     def __init__(self, t_ax, dta_avg) -> None:
#         self.t_ax = t_ax
#         self.dta_avg = dta_avg 

#     def transient_plot2D(self):
#         fig, ax = plt.subplots(1, 1)
#         ax.pcolormesh(np.real(self.dta_avg).T)
#         fig.show()
    
#     def plot(self):
#         fig, ax = plt.subplots(1, 1)
#         ax.plot(self.dta_x[0], np.real(self.dta_ev), label='Re')
#         ax.plot(self.dta_x[0], np.imag(self.dta_ev), label='Im')
#         ax.legend()
#         fig.show()
