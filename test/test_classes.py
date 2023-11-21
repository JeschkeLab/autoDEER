import numpy as np
from autodeer.classes import *
import pytest
from autodeer.pulses import RectPulse

# =============================================================================
# Interface
# =============================================================================

def test_make_interface():
    Interface1 = Interface()


# =============================================================================
# Parameters
# =============================================================================

def test_make_parameter():
    Par1 = Parameter(
            name="Par1", value=10, unit="us", description="The first parameter")

@pytest.fixture
def Par1():
    Par1 = Parameter(
            name="Par1", value=10, unit="us", description="The first parameter")
    return Par1

@pytest.fixture
def Par2():
    Par2 = Parameter(
        name="Par2", value=10,dim=5,step=2, unit="us", description="The second parameter")
    return Par2

def test_parameter_copy(Par1):
    Par2 = Par1.copy()
    assert Par2.value == 10
    assert Par2.unit == "us"
    assert Par2.description == "The first parameter"

def test_parameter_eq(Par1):
    Par2 = Par1.copy()
    assert Par1 == Par2

def test_parameter_add_scalar_number(Par1):
    Par2 = Par1 + 2
    assert Par2.value == 12


def test_parameter_sub_scalar_number(Par1):
    Par2 = Par1 - 2
    assert Par2.value == 8


def test_parameter_mul_scalar_number(Par1):
    Par2 = Par1 * 2
    assert Par2.value == 20

def test_parameter_add_static_param(Par1):
    Par3 = Parameter(
        name="Par3", value=2, unit="us", description="The first parameter")
    Par2 = Par1 + Par3
    assert Par2.value == 12


def test_parameter_sub_static_param(Par1):
    Par3 = Parameter(
        name="Par3", value=2, unit="us", description="The first parameter")
    Par2 = Par1 - Par3
    assert Par2.value == 8


def test_parameter_mul_static_param(Par1):
    Par3 = Parameter(
        name="Par3", value=2, unit="us", description="The first parameter")
    Par2 = Par1 * Par3
    assert Par2.value == 20

def test_parameter_add_dynamic_param(Par1):
    Par3 = Parameter(
        name="Par3", value=2, dim=5, step=1, unit="us", description="The first parameter")
    Par2 = Par1 + Par3
    assert Par2.value == 12
    assert np.array_equal(Par2.axis[0], Par3.axis[0])


def test_parameter_sub_dynamic_param(Par1):
    Par3 = Parameter(
        name="Par3", value=2, dim=5, step=1, unit="us", description="The first parameter")
    Par2 = Par1 - Par3
    assert Par2.value == 8
    assert np.array_equal(Par2.axis[0], Par3.axis[0])


def test_parameter_mul_dynamic_param(Par2):
    Par3 = Parameter(
        name="Par3", value=2, dim=5, step=1, unit="us", description="The first parameter")
    with pytest.raises(Exception) as e_info:
        Par = Par2 * Par3

def test_parameter_static_is_static(Par1):
    assert Par1.is_static() == True

def test_parameter_dynamic_is_static(Par2):
    assert Par2.is_static() == False

def test_parameter_remove_dynamic(Par2):
    Par2.remove_dynamic()
    assert Par2.is_static() == True

def test_parameter_get_axis(Par2):
    axis = Par2.get_axis()
    assert np.array_equal(axis, 10+np.arange(0,10,2))

def test_parameter_add_axis(Par1):
    Par1.add_axis(10,np.arange(0,10,2))
    assert np.array_equal(Par1.get_axis(), 10+np.arange(0,10,2))



    