import numpy as np
from autodeer.classes import *
import pytest
from autodeer.pulses import RectPulse


# =============================================================================
# Parameters
# =============================================================================

def test_make_parameter():
    Par1 = Parameter(
            name="Par1", value=10, unit="us", description="The first parameter")

Par1 = Parameter(
            name="Par1", value=10, unit="us", description="The first parameter")
Par1.add_axis(0,[0,2,4,6,8,10])

def test_parameter_add_scalar_number():
    Par2 = Par1 + 2
    assert Par2.value == 12
    assert np.array_equal(Par2.axis[0], [0,2,4,6,8,10])


def test_parameter_sub_scalar_number():
    Par2 = Par1 - 2
    assert Par2.value == 8
    assert np.array_equal(Par2.axis[0], [0,2,4,6,8,10])


def test_parameter_mul_scalar_number():
    Par2 = Par1 * 2
    assert Par2.value == 20
    assert np.array_equal(Par2.axis[0], [0,4,8,12,16,20])

def test_parameter_add_static_param():
    Par3 = Parameter(
        name="Par3", value=2, unit="us", description="The first parameter")
    Par2 = Par1 + Par3
    assert Par2.value == 12
    assert np.array_equal(Par2.axis[0], [0,2,4,6,8,10])


def test_parameter_sub_static_param():
    Par3 = Parameter(
        name="Par3", value=2, unit="us", description="The first parameter")
    Par2 = Par1 - Par3
    assert Par2.value == 8
    assert np.array_equal(Par2.axis[0], [0,2,4,6,8,10])


def test_parameter_mul_static_param():
    Par3 = Parameter(
        name="Par3", value=2, unit="us", description="The first parameter")
    Par2 = Par1 * Par3
    assert Par2.value == 20
    assert np.array_equal(Par2.axis[0], [0,4,8,12,16,20])

def test_parameter_add_dynamic_param():
    Par3 = Parameter(
        name="Par3", value=2, unit="us", description="The first parameter")
    Par3.add_axis(0,[0,1,2,3,4,5])
    Par2 = Par1 + Par3
    assert Par2.value == 12
    assert np.array_equal(Par2.axis[0], [0, 3, 6, 9, 12, 15])


def test_parameter_sub_dynamic_param():
    Par3 = Parameter(
        name="Par3", value=2, unit="us", description="The first parameter")
    Par3.add_axis(0,[0,1,2,3,4,5])
    Par2 = Par1 - Par3
    assert Par2.value == 8
    assert np.array_equal(Par2.axis[0], [0, 1, 2, 3, 4, 5])


def test_parameter_mul_dynamic_param():
    Par3 = Parameter(
        name="Par3", value=2, unit="us", description="The first parameter")
    Par3.add_axis(0,[0,1,2,3,4,5])
    with pytest.raises(Exception) as e_info:
        Par2 = Par1 * Par3


# =============================================================================
# Sequences
# =============================================================================


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

    