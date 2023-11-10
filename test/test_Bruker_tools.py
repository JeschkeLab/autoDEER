from autodeer.hardware.Bruker_tools import *
from autodeer.hardware.Bruker_MPFU import _MPFU_channels

from autodeer.sequences import *
import numpy as np
import pytest

def test_PulseSpel_MissingScale():
    B = 12200
    LO = 34.0
    reptime = 3e3
    averages = 1
    shots = 100
    test_sequence = HahnEchoSequence(B=B, LO=LO, reptime=reptime, 
        averages=averages, shots=shots
    )
    test_sequence.addPulsesProg(
        [1,2],
        ["t","t"],
        0,
        np.arange(200,5000,100)
    )
    test_sequence.convert()
    with pytest.raises(RuntimeError):   
        pulse_spel = PulseSpel(test_sequence, AWG=True)

def test_PulseSpel_HahnEcho():
    B = 12200
    LO = 34.0
    reptime = 3e3
    averages = 1
    shots = 100
    test_sequence = HahnEchoSequence(B=B, LO=LO, reptime=reptime, 
        averages=averages, shots=shots
    )
    test_sequence.addPulsesProg(
        [1,2],
        ["t","t"],
        0,
        np.arange(200,5000,100)
    )
    test_sequence.pulses[0].scale.value = 0.1
    test_sequence.pulses[1].scale.value = 0.1
    test_sequence.convert()
    pulse_spel = PulseSpel(test_sequence, AWG=True)
    _MPFU_channels(test_sequence)
    pulse_spel = PulseSpel(test_sequence, MPFU=["+<x>","-<x>","+<y>","-<y>"], AWG=False)

def test_PulseSpel_DEER():
    tau1 = 400
    tau2 = 2000
    dt = 16
    B = 12200
    LO = 34.0
    reptime = 3e3
    averages = 1
    shots = 100
    tp = 12
    test_sequence = DEERSequence(
        tau1=tau1, tau2=tau2, dt=dt, B=B, LO=LO, reptime=reptime, 
        averages=averages, shots=shots)

    test_sequence.four_pulse(tp)
    test_sequence.pulses[0].scale.value = 0.1
    test_sequence.pulses[1].scale.value = 0.1
    test_sequence.pulses[2].scale.value = 0.1
    test_sequence.pulses[3].scale.value = 0.1

    test_sequence.convert()
    pulse_spel = PulseSpel(test_sequence, AWG=True)
    _MPFU_channels(test_sequence)
    pulse_spel = PulseSpel(test_sequence, MPFU=["+<x>","-<x>","+<y>","-<y>"], AWG=False)

@pytest.fixture()
def build_deer_sequence():
    deer = Sequence(name="4p DEER", B=12220, LO=34.0,reptime=3e3, averages=1, shots=100)

    # tau1 = Parameter(name="tau1", value=400, unit="ns", step=16, dim=8, description="The first interpulse delays")
    tau1 = Parameter(name="tau1", value=400, unit="ns", description="The first interpulse delays")

    tau2 = Parameter(name="tau2", value=2500, unit="ns", description="The second interpulse delays")
    t = Parameter(name="t", value=-160, step=24, dim=120, unit="ns", description="The time axis")

    exc_pulse = RectPulse(freq=0, tp=12, flipangle=np.pi/2)
    ref_pulse = RectPulse(freq=0, tp=12, flipangle=np.pi)
    pump_pulse = RectPulse(freq=-0.1, tp=12, flipangle=np.pi)

    deer.addPulse(exc_pulse.copy(t=0, pcyc={"phases":[0, np.pi], "dets":[1,-1]}))
    deer.addPulse(ref_pulse.copy(t=tau1))
    deer.addPulse(pump_pulse.copy(t=2*tau1+t))
    deer.addPulse(ref_pulse.copy(t=2*tau1+tau2))
    deer.addPulse(Detection(t=2*(tau1+tau2), tp=512))

    deer.evolution([t,tau1], reduce=[tau1])
    return deer


def test_write_pulsespel_file_deer_MPFU(build_deer_sequence):
    seq = build_deer_sequence
    channels = _MPFU_channels(seq)
    MPFU = ["+<x>","-<x>","+<y>","-<y>"]
    def_file, exp_file = write_pulsespel_file(seq,AWG=False,MPFU=MPFU)
    print(exp_file)

def test_write_pulsespel_file_respro_MPFU():
    seq = ResonatorProfileSequence(
        B=12200, LO=34.0, reptime=3e3, averages=1, shots=10)
    channels = _MPFU_channels(seq)
    MPFU = ["+<x>","-<x>","+<y>","-<y>"]
    def_file, exp_file = write_pulsespel_file(seq,AWG=False,MPFU=MPFU)
    print(exp_file)
    




def test_PulseSpel_ResonatorProfile():

    pass
def test_PulseSpel_TWTProfile():
    pass