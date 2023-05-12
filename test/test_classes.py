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

def test_parameter_to_from_json():
    json_file = Par1._to_json()
    new_par = Par1._from_json(json_file)
    print(new_par)


# =============================================================================
# Sequences
# =============================================================================



    