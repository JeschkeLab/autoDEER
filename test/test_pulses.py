from autodeer.pulses import *
import pytest

def test_pulse_init():
    pulse = Pulse(tp=10, scale=0.5)
    assert pulse.tp.value == 10
    assert pulse.scale.value == 0.5

def test_pulse_add_phase_cycle():
    pulse = Pulse(tp=10, scale=0.5)
    pulse._addPhaseCycle([0, 90, 180, 270])
    assert pulse.pcyc["Phases"] == [0, 90, 180, 270]
    assert pulse.pcyc["DetSigns"] == [1, 1, 1, 1]


def test_pulse_is_static():
    pulse = Pulse(tp=10, scale=0.5)
    assert pulse.is_static()
    pulse = Pulse(tp=10, scale=0.5, t=0)
    assert pulse.is_static()
    pulse = Pulse(tp=10, scale=0.5, t=0, pcyc=None)
    assert pulse.is_static()
    tp = Parameter(
                name="tp", value=10, unit="us", description="The pulse length",
                axis=np.arange(0,10,1), axis_id=0)
    pulse = Pulse(tp=tp, scale=0.5, t=0, pcyc=None)
    assert pulse.is_static()

def test_pulse_is_delay_focused():
    pulse = Pulse(tp=10, scale=0.5)
    assert pulse.isDelayFocused()
    pulse = Pulse(tp=10, scale=0.5, t=0)
    assert not pulse.isDelayFocused()

def test_pulse_is_pulse_focused():
    pulse = Pulse(tp=10, scale=0.5)
    assert not pulse.isPulseFocused()
    pulse = Pulse(tp=10, scale=0.5, t=0)
    assert pulse.isPulseFocused()

def test_pulse_plot():
    pulse = Pulse(tp=10, scale=0.5)
    pulse._buildFMAM(lambda x: (np.ones_like(x), np.zeros_like(x)))
    pulse.plot()

def test_pulse_exciteprofile():
    pulse = Pulse(tp=10, scale=0.5, flipangle=np.pi/2)
    pulse.func = lambda x: (np.ones_like(x), np.zeros_like(x))
    freqs = np.linspace(-.25, .25, 1001)
    Mx,My,Mz = pulse.exciteprofile(freqs=freqs)
    assert len(Mx) == 1001
    assert len(My) == 1001
    assert len(Mz) == 1001

def test_rectpulse_init():
    pulse = RectPulse(tp=10, freq=1, t=0, flipangle=np.pi/2)
    assert pulse.tp.value == 10
    assert pulse.freq.value == 1
    assert pulse.t.value == 0
    assert pulse.flipangle.value == np.pi/2

def test_rectpulse_func():
    pulse = RectPulse(tp=10, freq=1, t=0, flipangle=np.pi/2)
    AM, FM = pulse.func(np.arange(10))
    assert np.allclose(AM, np.ones(10))
    assert np.allclose(FM, np.ones(10) )
    
def test_hspulse_init():
    pulse = HSPulse(tp=10, order1=1, order2=6, beta=20, BW=1, init_freq=0)
    assert pulse.tp.value == 10
    assert pulse.order1.value == 1
    assert pulse.order2.value == 6
    assert pulse.beta.value == 20
    assert pulse.BW.value == 1
    assert pulse.init_freq.value == 0
    assert pulse.final_freq.value == 1

def test_hspulse_func():
    pulse = HSPulse(tp=128, order1=1, order2=6, beta=20, BW=0.1, init_freq=0)
    AM, FM = pulse.func(np.arange(10))
    # assert np.allclose(AM, np.ones(10))
    assert np.allclose(FM.min(), 0)
    assert np.allclose(FM.max(), 0.1)

