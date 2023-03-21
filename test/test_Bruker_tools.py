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


def test_PulseSpel_ResonatorProfile():

    pass
def test_PulseSpel_TWTProfile():
    pass