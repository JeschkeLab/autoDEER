from autodeer.sequences import *
import pytest
from matplotlib.figure import Figure

# ------------------ Test the Sequence class ------------------

def test_init_sequence():
    seq = Sequence(name="name", B=12220, LO=34.0,reptime=3e3, averages=1, shots=100)
    assert seq.name == "name"
    assert seq.B.value == 12220
    assert seq.LO.value == 34.0
    assert seq.reptime.value == 3e3
    assert seq.averages.value == 1
    assert seq.shots.value == 100

def test_add_pulse():
    seq = Sequence(name="name", B=12220, LO=34.0,reptime=3e3, averages=1, shots=100)
    pulse = RectPulse(t=0, freq=0, tp=12, flipangle=np.pi/2)
    seq.addPulse(pulse)
    assert seq.pulses[0] == pulse
    assert len(seq.pulses) == 1

def test_add_pulse_list():
    seq = Sequence(name="name", B=12220, LO=34.0,reptime=3e3, averages=1, shots=100)
    pulse = RectPulse(t=0,freq=0, tp=12, flipangle=np.pi/2)
    seq.addPulse([pulse, pulse])
    assert seq.pulses[0] == pulse
    assert seq.pulses[1] == pulse
    assert len(seq.pulses) == 2
    
def test_add_pulse_evolution():
    seq=Sequence(name="name", B=12220, LO=34.0,reptime=3e3, averages=1, shots=100)
    t = Parameter(name="t", value=-160, step=24, dim=120, unit="ns", description="The time axis")
    pulse = RectPulse(t=t, freq=0, tp=12, flipangle=np.pi/2)
    seq.addPulse(pulse)
    seq.evolution([t])
    assert seq.progTable["Variable"][0] == "t"

def test_adjust_step():
    seq=Sequence(name="name", B=12220, LO=34.0,reptime=3e3, averages=1, shots=100)
    t = Parameter(name="t", value=-161, step=25, dim=120, unit="ns", description="The time axis")
    pulse = RectPulse(t=t, freq=0, tp=12, flipangle=np.pi/2)
    seq.addPulse(pulse)
    seq.evolution([t])
    seq.adjust_step(2)

    assert np.equal(np.diff(seq.progTable["axis"][0]),24).all()
    assert np.equal(seq.pulses[0].t.value, -160)

def test_adjust_step_linked_params():
    # check that adujst step also works with linked parameters
    seq=Sequence(name="name", B=12220, LO=34.0,reptime=3e3, averages=1, shots=100)
    t = Parameter(name="t", value=151, step=25, dim=120, unit="ns", description="The time axis")
    pulse1 = RectPulse(t=t, freq=0.5, tp=12, flipangle=np.pi/2)
    pulse2 = RectPulse(t=2*t, freq=0.5, tp=12, flipangle=np.pi/2)
    seq.addPulse(pulse1)
    seq.addPulse(pulse2)
    seq.evolution([t])
    seq.adjust_step(2)

    assert np.equal(np.diff(seq.progTable["axis"][0]),24).all()
    assert np.equal(np.diff(seq.progTable["axis"][1]),48).all()
    assert np.equal(seq.pulses[0].t.value, 152)
    assert np.equal(seq.pulses[1].t.value, 304)
    assert np.equal(seq.pulses[1].freq.value, 0.5)


@pytest.fixture
def test_sequence():
    seq=Sequence(name="name", B=12220, LO=34.0,reptime=3e3, averages=1, shots=100)
    t = Parameter(name="t", value=-160, step=24, dim=120, unit="ns", description="The time axis")
    pulse = RectPulse(t=t, freq=0, tp=12, flipangle=np.pi/2,pcyc={"phases":[0, np.pi], "dets":[1,-1]})
    seq.addPulse(pulse)
    seq.evolution([t])
    return seq

def test_estimate_time(test_sequence):
    assert test_sequence._estimate_time() == pytest.approx(72, abs=0.5)

def test_print(test_sequence):
    assert isinstance(test_sequence.__str__(),str)

def test_plot_pulse_exc(test_sequence):
    fig = test_sequence.plot_pulse_exc()
    assert isinstance(fig, Figure)

def test_to_dict(test_sequence):
    dic = test_sequence._to_dict()
    assert isinstance(test_sequence._to_dict(), dict)
    assert len(dic.pulses) == len(test_sequence.pulses)

def test_sequence_to_from_json(test_sequence):
    json_file = test_sequence._to_json()
    new_seq = Sequence._from_json(json_file)

    assert new_seq.name == test_sequence.name
    assert new_seq.B.value == test_sequence.B.value
    assert new_seq.LO.value == test_sequence.LO.value
    assert new_seq.reptime.value == test_sequence.reptime.value
    assert new_seq.averages.value == test_sequence.averages.value
    assert new_seq.shots.value == test_sequence.shots.value
    assert new_seq.pulses[0].t.value == test_sequence.pulses[0].t.value
    assert new_seq.pulses[0].tp.value == test_sequence.pulses[0].tp.value


# ---------------Test the pre-built sequences ------------------

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
    tau1 = 0.4
    tau2 = 2
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
    test_sequence.select_pcyc('DC')
    
    assert test_sequence.pulses[0].t.value == 0
    assert test_sequence.pulses[0].tp.value == tp
    assert test_sequence.pulses[1].t.value == tau1*1e3
    assert test_sequence.pulses[1].tp.value == tp
    assert test_sequence.pulses[2].tp.value == tp
    assert test_sequence.pulses[3].t.value == (2*tau1 + tau2)*1e3
    assert test_sequence.pulses[3].tp.value == tp
    assert test_sequence.pulses[4].t.value == (2*tau1 + 2*tau2)*1e3

    assert test_sequence.pulses[0].pcyc == {'Phases': [0,np.pi], 'DetSigns': [1, -1]}

    assert test_sequence.B.value == B
    assert test_sequence.LO.value == LO
    assert test_sequence.reptime.value == reptime
    assert test_sequence.averages.value == averages
    assert test_sequence.shots.value == shots

def test_DEER_5p():
    tau1 = 2
    tau2 = tau1
    tau3 = 0.3
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
    test_sequence.select_pcyc('16step_5p')

    assert test_sequence.pulses[0].t.value == 0
    assert test_sequence.pulses[0].tp.value == tp
    assert test_sequence.pulses[1].t.value == (tau1-tau3)*1e3
    assert test_sequence.pulses[1].tp.value == tp
    assert test_sequence.pulses[2].t.value == tau1*1e3
    assert test_sequence.pulses[2].tp.value == tp
    assert test_sequence.pulses[3].tp.value == tp
    assert test_sequence.pulses[4].t.value == (2*tau1 + tau2)*1e3
    assert test_sequence.pulses[4].tp.value == tp
    assert test_sequence.pulses[5].t.value == (2*tau1 + 2*tau2)*1e3

    assert test_sequence.pulses[2].pcyc == {'Phases': [0, np.pi/2, np.pi, -np.pi/2], 'DetSigns': [1, -1, 1, -1]}
    assert test_sequence.pulses[3].pcyc == {'Phases': [0, np.pi/2, np.pi, -np.pi/2], 'DetSigns': [1, 1, 1, 1]}

    assert test_sequence.B.value == B
    assert test_sequence.LO.value == LO
    assert test_sequence.reptime.value == reptime
    assert test_sequence.averages.value == averages
    assert test_sequence.shots.value == shots

def test_DEER_7p():
    tau1 = 2
    tau2 = tau1
    tau3 = 0.3
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

    test_sequence.seven_pulse(tp)
    test_sequence.select_pcyc('32step_7p')

    assert len(test_sequence.pulses) == 8

    assert test_sequence.pulses[1].pcyc == {'Phases': [0, np.pi/2, np.pi, -np.pi/2], 'DetSigns': [1, -1, 1, -1]}
    assert test_sequence.pulses[3].pcyc == {'Phases': [0, np.pi,], 'DetSigns': [1, 1]}
    assert test_sequence.pulses[4].pcyc == {'Phases': [0, np.pi,], 'DetSigns': [1, 1]}
    assert test_sequence.pulses[5].pcyc == {'Phases': [0, np.pi,], 'DetSigns': [1, 1]}


def test_nDEER_CP():
    tau1 = 2
    tau2 = tau1
    tau3 = 0.3
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

    test_sequence.nDEER_CP(4)

    assert len(test_sequence.pulses) == 7


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
    assert test_sequence.pulses[0].tp.value == 0
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

def test_ReptimeScan():
    B = 12200
    LO = 34.0
    reptime = 3e4
    averages = 1
    shots = 100
    test_sequence = ReptimeScan(B=B, LO=LO,reptime=3e3, reptime_max=reptime, 
        averages=averages, shots=shots
    )

    assert test_sequence.progTable["Variable"][0] == "reptime"
    axis = test_sequence.reptime.get_axis()
    assert axis.min() == 0
    assert axis.max() == pytest.approx(reptime, abs=1e3)
    assert axis.shape[0] == 100
    assert test_sequence.reptime.value == 3e3


def test_CarrPurcellSequence():
    B = 12200
    LO = 34.0
    reptime = 3e3
    averages = 1
    shots = 100
    test_sequence = CarrPurcellSequence(B=B, LO=LO, reptime=reptime, 
        averages=averages, shots=shots,n=2, tau=10
    )

    assert len(test_sequence.progTable["Variable"]) == 3
    assert test_sequence.progTable["Variable"] == ['t','t','t']
