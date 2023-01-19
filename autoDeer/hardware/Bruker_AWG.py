from autodeer.openepr import Interface, dataset, Sequence, Delay, Detection
from autodeer.hardware.xepr_api_adv import xepr_api
from autodeer.hardware.Bruker_tools import PulseSpel, run_general
from autodeer.sequences import HahnEchoSequence

import tempfile
import time
from scipy.optimize import minimize_scalar, curve_fit
import numpy as np

# =============================================================================


class BrukerAWG(Interface):

    def __init__(self, config_file:str, d0=600) -> None:

        self.api = xepr_api(config_file)
        self.spec_config = self.api.config["Spectrometer"]
        self.bridge_config = self.api.spec_config["Bridge"]
        
        self.temp_dir = tempfile.mkdtemp("autoDEER")

        self.d0 = d0
        
        super().__init__()


    def connect(self) -> None:

        self.api.connect()
        return super().connect()

    def acquire_dataset(self) -> dataset:
        return self.api.acquire_dataset()
    
    def launch(self, sequence, savename: str, start=True):
        
        PSpel = PulseSpel(sequence, AWG=True)

        PSpel_file = self.temp_dir + "/autoDEER_PulseSpel"
        PSpel.save(PSpel_file)

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
    
    def tune(self, sequence) -> None:

        self.api.set_field(sequence.B.value)
        self.api.set_freq(sequence.LO.value)

        pass
    
    def tune_nutation(self, test_pulse):
        Nutation = open.Sequence(
            name="Nutation Tuning", B=12248.7, LO=33.970715,
            reptime=3000, averages=1, shots=100)

        pi2_pulse = open.RectPulse(t=0, tp=12, freq=0, pcyc = [0, np.pi], flipangle=np.pi/2)
        pi_pulse = open.RectPulse(t=500, tp=12, freq=0, pcyc= [0], flipangle=np.pi)

        Nutation.addPulse(test_pulse)
        Nutation.addPulse(pi2_pulse)
        Nutation.addPulse(pi_pulse)

        Nutation.addPulse(open.Detection(t=1e3, tp=48))
        Nutation.convert()
        pass

    
    def isrunning(self) -> bool:
        return self.api.is_exp_running()

    def terminate(self) -> None:
        return self.api.abort_exp()