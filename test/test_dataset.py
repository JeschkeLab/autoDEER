from autodeer.classes import *
from autodeer.dataset import *
from autodeer.sequences import *
import numpy as np
from scipy.signal import hilbert

def test_make_dataset():
    in_data = np.random.random(100)
    in_axes = [np.linspace(0,100,100)]
    test_dataset = Dataset(axes=in_axes, data=in_data, force_complex=False)

    assert np.array_equal(test_dataset.data, in_data)
    assert np.array_equal(test_dataset.axes, in_axes)

def test_hilbert():
    in_data = np.random.random(100)
    in_axes = [np.linspace(0,100,100)]
    test_dataset = Dataset(axes=in_axes, data=in_data, force_complex=True)

    assert np.array_equal(test_dataset.data, hilbert(in_data))

def test_save():
    in_data = np.random.random(100)
    in_axes = [np.linspace(0,100,100)]
    test_dataset = Dataset(axes=in_axes, data=in_data, force_complex=True)
    test_dataset.save("test_dataset.json")

    new_dataset = Dataset.load("test_dataset.json")
    assert np.array_equal(test_dataset.data, new_dataset.data)
    assert np.array_equal(test_dataset.axes, new_dataset.axes)
