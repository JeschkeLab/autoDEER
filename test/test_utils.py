from autodeer.utils import *

def test_gcd():
    assert gcd([2.0, 4.0, 6.0, 8.0]) == 2
    assert gcd([1.5, 3.0]) == 1.5
    assert gcd([3, 6, 9, 12]) == 3
    assert gcd([5, 10, 15, 20]) == 5
    assert gcd([7, 14, 21, 28]) == 7
    assert gcd([8, 12, 16, 20]) == 4
    assert gcd([10, 20, 30, 40]) == 10
    assert gcd([15, 25, 35, 45]) == 5
    assert gcd([18, 24, 30, 36]) == 6
    assert gcd([21, 28, 35, 42]) == 7
    assert gcd([27, 36, 45, 54]) == 9

def test_transpose_dict_of_list():
    d = {'a': [1, 2, 3], 'b': [4, 5, 6], 'c': [7, 8, 9]}
    expected = [{'a': 1, 'b': 4, 'c': 7}, {'a': 2, 'b': 5, 'c': 8}, {'a': 3, 'b': 6, 'c': 9}]
    assert transpose_dict_of_list(d) == expected

    d = {'a': [1, 2], 'b': [3, 4], 'c': [5, 6]}
    expected = [{'a': 1, 'b': 3, 'c': 5}, {'a': 2, 'b': 4, 'c': 6}]
    assert transpose_dict_of_list(d) == expected

    d = {'a': [1], 'b': [2], 'c': [3]}
    expected = [{'a': 1, 'b': 2, 'c': 3}]
    assert transpose_dict_of_list(d) == expected

    d = {'a': [], 'b': [], 'c': []}
    expected = []
    assert transpose_dict_of_list(d) == expected

def test_transpose_list_of_dicts():
    d = [{'a': 1, 'b': 4, 'c': 7}, {'a': 2, 'b': 5, 'c': 8}, {'a': 3, 'b': 6, 'c': 9}]
    expected = {'a': [1, 2, 3], 'b': [4, 5, 6], 'c': [7, 8, 9]}
    assert transpose_list_of_dicts(d) == expected

    d = [{'a': 1, 'b': 3, 'c': 5}, {'a': 2, 'b': 4, 'c': 6}]
    expected = {'a': [1, 2], 'b': [3, 4], 'c': [5, 6]}
    assert transpose_list_of_dicts(d) == expected

    d = [{'a': 1, 'b': 2, 'c': 3}]
    expected = {'a': [1], 'b': [2], 'c': [3]}
    assert transpose_list_of_dicts(d) == expected

    d = []
    expected = {}
    assert transpose_list_of_dicts(d) == expected

import numpy as np
from autodeer import Parameter
import pytest

def test_val_in_ns():
    # Test with no axis
    p = Parameter(name='p', value=10, unit="ns")
    assert val_in_ns(p) == 10
    p = Parameter(name='p', value=10, unit="us")
    assert val_in_ns(p) == 10000

    # Test with one axis
    p = Parameter(name='p', value=10, unit="ns", step=5, dim=1)
    assert np.all(val_in_ns(p) == np.array([10]))

    p = Parameter(name='p', value=10, unit="ns", step=5, dim=3)
    assert np.all(val_in_ns(p) == np.array([10,15,20]))

    # Test with two axis
    p = Parameter(name='p', value=10, unit="ns", step=5, dim=2)
    p.add_axis(2,axis=np.array([1,2]))
    
    with pytest.raises(ValueError):
        val_in_ns(p)

def test_val_in_us():
    # Test with no axis
    p = Parameter(name='p', value=10, unit="ns")
    assert val_in_us(p) == 0.01
    p = Parameter(name='p', value=10, unit="us")
    assert val_in_us(p) == 10

    # Test with one axis
    p = Parameter(name='p', value=10, unit="us", step=5, dim=1)
    assert np.all(val_in_us(p) == np.array([10]))

    p = Parameter(name='p', value=10, unit="us", step=5, dim=3)
    assert np.all(val_in_us(p) == np.array([10,15,20]))

    # Test with two axis
    p = Parameter(name='p', value=10, unit="us", step=5, dim=2)
    p.add_axis(2,axis=np.array([1,2]))
    
    with pytest.raises(ValueError):
        val_in_us(p)
