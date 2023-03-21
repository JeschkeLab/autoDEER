from autodeer.sequences import *
import pytest

def test_DEER_4p():
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
    
    assert test_sequence.pulses[0].t.value == 0
    assert test_sequence.pulses[0].tp.value == tp
    assert test_sequence.pulses[1].t.value == tau1
    assert test_sequence.pulses[1].tp.value == tp
    assert test_sequence.pulses[2].tp.value == tp
    assert test_sequence.pulses[3].t.value == 2*tau1 + tau2
    assert test_sequence.pulses[3].tp.value == tp
    assert test_sequence.pulses[4].t.value == 2*tau1 + 2*tau2

    assert test_sequence.B.value == B
    assert test_sequence.LO.value == LO
    assert test_sequence.reptime.value == reptime
    assert test_sequence.averages.value == averages
    assert test_sequence.shots.value == shots

def test_DEER_5p():
    tau1 = 2000
    tau2 = tau1
    tau3 = 200
    dt = 16
    B = 12200
    LO = 34.0
    reptime = 3e3
    averages = 1
    shots = 100
    tp = 12
    test_sequence = DEERSequence(
        tau1=tau1, tau2=tau2, tau3=tau3, dt=dt, B=B, LO=LO, reptime=reptime, 
        averages=averages, shots=shots)

    test_sequence.five_pulse(tp)
    
    assert test_sequence.pulses[0].t.value == 0
    assert test_sequence.pulses[0].tp.value == tp
    assert test_sequence.pulses[2].t.value == tau1
    assert test_sequence.pulses[2].tp.value == tp
    assert test_sequence.pulses[3].tp.value == tp
    assert test_sequence.pulses[4].t.value == 2*tau1 + tau2
    assert test_sequence.pulses[4].tp.value == tp
    assert test_sequence.pulses[5].t.value == 2*tau1 + 2*tau2

    assert test_sequence.B.value == B
    assert test_sequence.LO.value == LO
    assert test_sequence.reptime.value == reptime
    assert test_sequence.averages.value == averages
    assert test_sequence.shots.value == shots

def test_HahnEchoSequence():
    B = 12200
    LO = 34.0
    reptime = 3e3
    averages = 1
    shots = 100
    test_sequence = HahnEchoSequence(B=B, LO=LO, reptime=reptime, 
        averages=averages, shots=shots
    )
    assert test_sequence.pulses[0].t.value == 0
    assert test_sequence.pulses[0].tp.value == 12
    assert test_sequence.pulses[1].t.value == 500
    assert test_sequence.pulses[1].tp.value == 12
    assert test_sequence.pulses[2].t.value == 1000

def test_ResonatorProfile():
    B = 12200
    LO = 34.0
    reptime = 3e3
    averages = 1
    shots = 100
    test_sequence = ResonatorProfileSequence(B=B, LO=LO, reptime=reptime, 
        averages=averages, shots=shots
    )

    assert test_sequence.pulses[0].t.value == 0
    assert test_sequence.pulses[0].tp.value == 4
    assert test_sequence.pulses[1].t.value == 2000
    assert test_sequence.pulses[1].tp.value == 16
    assert test_sequence.pulses[2].t.value == 2500
    assert test_sequence.pulses[2].tp.value == 32
    assert test_sequence.pulses[3].t.value == 3000

    assert test_sequence.progTable["Variable"][0] == "tp"
    assert test_sequence.progTable["axis"][0].min() == 0
    assert test_sequence.progTable["axis"][0].max() == 64
    
def test_TWTProfile():
    B = 12200
    LO = 34.0
    reptime = 3e3
    averages = 1
    shots = 100
    test_sequence = TWTProfileSequence(B=B, LO=LO, reptime=reptime, 
        averages=averages, shots=shots
    )

    assert test_sequence.pulses[0].t.value == 0
    assert test_sequence.pulses[0].tp.value == 4
    assert test_sequence.pulses[1].t.value == 2000
    assert test_sequence.pulses[1].tp.value == 16
    assert test_sequence.pulses[2].t.value == 2500
    assert test_sequence.pulses[2].tp.value == 32
    assert test_sequence.pulses[3].t.value == 3000

    assert test_sequence.progTable["Variable"][0] == "tp"
    assert test_sequence.progTable["Variable"][1] == "scale"
    assert test_sequence.progTable["axis"][0].min() == 0
    assert test_sequence.progTable["axis"][0].max() == 64
    assert test_sequence.progTable["axis"][1].min() == pytest.approx(0)
    assert test_sequence.progTable["axis"][1].max() == pytest.approx(1)
