from autodeer.classes import *
from autodeer.dataset import *
from autodeer.sequences import *
import numpy as np
from scipy.signal import hilbert
import pytest

def test_create_dataset_from_sequence():
    seq = FieldSweepSequence(12200, 500, 34.4, reptime=3e3, averages=1, shots=100)
    data = np.ones(seq.B.axis[0]['axis'].shape[0])
    dset = create_dataset_from_sequence(data,seq)
    
def test_create_dataset_from_axis():
    t = np.linspace(0,3,50)
    V = np.ones(t.shape[0])
    params = {"par1":1, "par2":2}
    dset = create_dataset_from_axes(V, t, params, axes_labels=['t'])

    assert np.array_equal(dset.data, V)
    assert np.array_equal(dset.t, t)
    assert dset.par1 == 1
    assert dset.par2 == 2
    

def test_create_dataset_from_bruker():
    file_path = "test/test_data/test_DEER.DTA"
    dset = create_dataset_from_bruker(file_path)

    assert hasattr(dset, "t")
    assert dset.t.shape[0] == 307
    assert dset.shape[0] == 307
    assert dset.B == 12102.45
    assert dset.LO == 34.00708

@pytest.fixture 
def dataset():
    file_path = "test/test_data/test_DEER.DTA"
    dset = create_dataset_from_bruker(file_path)
    return dset

@pytest.mark.skip(reason="Not implemented")
def test_save_dataset(dataset):
    dataset.epr.save("test_dataset")
    

def test_correct_phase(dataset):
    dset = dataset.epr.correctphase
    assert not np.iscomplexobj(dset.data)

def test_correct_phase_full_outoput(dataset):
    dset = dataset.epr.correctphasefull
    assert np.iscomplexobj(dset.data)
    

def test_calc_snr(dataset):
    dset = dataset.epr.correctphase
    dset_snr = dset.epr.SNR
    assert np.abs(dset_snr - 96) / 96 < 0.1
    pass
