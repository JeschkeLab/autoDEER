from autodeer.sequences import *
import pytest


def test_build_DEER_sequence():
    deer = Sequence(name="4p DEER", B=12220, LO=34.0,reptime=3e3, averages=1, shots=100)

    tau1 = Parameter(name="tau1", value=400, unit="ns", step=16, dim=8, description="The first interpulse delays")
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

    print(deer)


def test_sequence_to_from_json():
    deer = Sequence(name="4p DEER", B=12220, LO=34.0,reptime=3e3, averages=1, shots=100)

    tau1 = Parameter(name="tau1", value=400, unit="ns", step=16, dim=8, description="The first interpulse delays")
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

    json_file = deer._to_json()
    new_deer = Sequence._from_json(json_file)
    print(new_deer)

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
