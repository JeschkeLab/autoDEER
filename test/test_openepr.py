import autodeer.openepr as openepr
import numpy as np


def test_Pulse():

    pulse = openepr.Pulse(0, 32, 0.5)
    assert pulse.t.value == 0
    assert pulse.tp.value == 32
    assert pulse.scale.value == 0.5


def test_RectPulse():

    pulse = openepr.RectPulse(0, 32, 0.1, 0.5)
    assert pulse.t.value == 0
    assert pulse.tp.value == 32
    assert pulse.scale.value == 0.5
    assert pulse.freq.value == 0.1

    assert np.array_equal(pulse.AM, np.ones(1000))
    assert np.array_equal(pulse.FM, np.zeros(1000))


def test_phasecycle():
    pulse = openepr.Pulse(0, 32, 0.5)
    pulse._addPhaseCycle([0, np.pi/2, np.pi, -np.pi/2])


def test_sequence():
    Hahn_echo = openepr.Sequence()
    Hahn_echo.addPulse(openepr.RectPulse(t=0, tp=16, freq=0, scale=0.1, flipangle=np.pi/2))
    Hahn_echo.addPulse(openepr.RectPulse(t=5e2, tp=32, freq=0.1, scale=0.1, flipangle=np.pi))
    Hahn_echo.addPulse(openepr.Detection(t=1e3,tp=512))

    Hahn_echo.addPulsesProg(
        pulses=[1,2],
        variables=['t','t'],
        axis_id=0,
        axis=np.arange(2e2,2e3,50),
    )

def build_hahn_sequence():
    Hahn_echo = openepr.Sequence()
    Hahn_echo.addPulse(openepr.RectPulse(t=0, tp=16, freq=0, scale=0.1, flipangle=np.pi/2))
    Hahn_echo.addPulse(openepr.RectPulse(t=5e2, tp=32, freq=0.1, scale=0.1, flipangle=np.pi))
    Hahn_echo.addPulse(openepr.Detection(t=1e3,tp=512))

    Hahn_echo.addPulsesProg(
        pulses=[1,2],
        variables=['t','t'],
        axis_id=0,
        axis=np.arange(2e2,2e3,50),
    )
def test_is_pulse_focused():
    Hahn_echo = build_hahn_sequence ()
    assert Hahn_echo.isPulseFocused() == True

def test_print_seq():
    Hahn_echo = build_hahn_sequence ()
    print(Hahn_echo)

def test_plot_pulse_exc():
    Hahn_echo = build_hahn_sequence ()
    Hahn_echo.plot_pulse_exc()
    
