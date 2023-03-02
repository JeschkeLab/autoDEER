import matlab.engine
from autodeer import Pulse, Delay, Detection, Sequence, RectPulse, ChirpPulse, HSPulse
from autodeer import HahnEchoSequence
import numpy as np
import os
import re
from autodeer import eprload
import time
from deerlab import correctphase


class ETH_awg_interface:
    """
    Represents the interface for connecting to Andrin Doll style spectrometers.
    """
    def __init__(self, awg_freq=1.5, dig_rate=2) -> None:
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
        
        self.connect()
        self.awg_freq = awg_freq
        self.dig_rate = dig_rate
        self.pulses = {}
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
                data = eprload(path)
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
        self.cur_exp = struct
        self.workspace['exp'] = self.cur_exp
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
        amp_tune.pulses[0].tp.value = tp
        amp_tune.pulses[0].scale.value = 0
        amp_tune.pulses[1].tp.value = tp * 2
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

            nut_tune.addPulsesProg(
                pulses=[0],
                variables=["scale"],
                axis_id = 0,
                axis= np.arange(0,0.9,0.02)
            )
            self.launch(nut_tune, "autoDEER_amptune", IFgain=1)

            while self.isrunning():
                time.sleep(10)
            dataset = self.acquire_dataset()
            data = correctphase(dataset.data)
            if data[np.abs(data).argmax()] < 0:
                data *= -1
            if np.isclose(pulse.flipangle.value, np.pi):
                scale = np.around(dataset.axes[0][data.argmin()],2)
            elif np.isclose(pulse.flipangle.value, np.pi/2):
                sign_changes = np.diff(np.sign(np.real(data)))
                scale = np.around(dataset.axes[0][np.nonzero(sign_changes)[0][0]],2)
            else:
                raise RuntimeError("Target pulse can only have a flip angle of either: ",
                                "pi or pi/2.")
            pulse.scale.value = scale
        
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

        struc["LO"] = round(float(sequence.LO.value - self.awg_freq), 3)
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
                vec = prog_table["axis"][i].astype(float)
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
                    vec = prog_table["axis"][i].astype(float)
                    centre_freq = (vec[-1] + vec[0])/2
                    LO = centre_freq - self.awg_freq
                    sequence.LO.value = LO
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

    def terminate_at(self, criterion, test_interval=10):
        """Terminates the experiment upon a specific condition being
        satisified. 

        Parameters
        ----------
        criterion : _type_
            The criteria to be tested.
        test_interval : int, optional
            How often should the criteria be tested in minutes, by default 10.
        """

        test_interval_seconds = test_interval * 60
        condition = False
        last_scan = 0

        start_time = time.time()

        while not condition:
            data = self.acquire_dataset()
            try:
                nAvgs = data.nAvgs.value
            except AttributeError:
                print("WARNING: Dataset missing number of averages(nAvgs)!")
                nAvgs = 1
            finally:
                if nAvgs < 1:
                    time.sleep(30)  # Replace with single scan time
                    continue
                elif nAvgs <= last_scan:
                    time.sleep(30)
                    continue

            condition = criterion.test(data)

            if not condition:
                end_time = time.time()

                if (end_time - start_time) < test_interval_seconds:
                    time.sleep(test_interval_seconds - (end_time - start_time))

        self.terminate()
        pass
