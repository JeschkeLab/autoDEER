import autoDeer.hardware.openepr as openepr
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
    Hahn_echo.addPulse(openepr.RectPulse(0, 16, 0, 0.1))
    Hahn_echo.addPulse(openepr.RectPulse(0, 32, 0, 0.1))
