from autodeer.classes import *
from autodeer.dataset import *
from autodeer.sequences import *
from autodeer.hardware.dummy import _simulate_field_sweep, _simulate_T2
import numpy as np
from scipy.signal import hilbert
import pytest

def test_create_dataset_from_sequence():
    seq = FieldSweepSequence(B=12200, Bwidth=500, LO=34.4, reptime=3e3, averages=1, shots=100)
    data = np.ones(seq.B.axis[0]['axis'].shape[0])
    extra_params = {"nAvgs":1}
    dset = create_dataset_from_sequence(data,seq, extra_params)

    assert 'autoDEER_Version' in dset.attrs
    assert 'nAvgs' in dset.attrs
    assert 'shots' in dset.attrs
    assert 'reptime' in dset.attrs
    assert 'nPcyc' in dset.attrs
    
def test_create_dataset_from_axis():
    t = np.linspace(0,3,50)
    V = np.ones(t.shape[0])
    params = {"par1":1, "par2":2}
    dset = create_dataset_from_axes(V, t, params, axes_labels=['t'])

    assert np.array_equal(dset.data, V)
    assert np.array_equal(dset.t, t)
    assert dset.par1 == 1
    assert dset.par2 == 2
    assert 'autoDEER_Version' in dset.attrs
    

def test_create_dataset_from_bruker():
    file_path = "test/test_data/test_DEER.DTA"
    dset = create_dataset_from_bruker(file_path)

    assert hasattr(dset, "t")
    assert dset.t.shape[0] == 307
    assert dset.shape[0] == 307
    assert dset.B == 12102.45
    assert dset.LO == 34.00708
    assert 'autoDEER_Version' in dset.attrs

@pytest.fixture 
def dataset_from_bruker():
    file_path = "test/test_data/test_DEER.DTA"
    dset = create_dataset_from_bruker(file_path)
    return dset

@pytest.fixture
def dataset_from_sequence():
    # Fix the seed for reproducibility
    rng = np.random.default_rng(seed=0)
    
    seq = FieldSweepSequence(B=12200, Bwidth=500, LO=34.4, reptime=3e3, averages=1, shots=100)

    t, data = _simulate_field_sweep(seq)
    data += rng.normal(0,0.04,data.shape) + 1j*rng.normal(0,0.04,data.shape)
    extra_params = {"nAvgs":1}
    dset = create_dataset_from_sequence(data,seq, extra_params)
    return dset
    

def test_correct_phase(dataset_from_sequence):
    dataset = dataset_from_sequence
    dset = dataset.epr.correctphase
    assert not np.iscomplexobj(dset.data)

def test_correct_phase_full_outoput(dataset_from_sequence):
    dataset = dataset_from_sequence
    dset = dataset.epr.correctphasefull
    assert np.iscomplexobj(dset.data)
    

def test_calc_snr(dataset_from_sequence):
    dataset = dataset_from_sequence
    dset = dataset.epr.correctphase
    dset_snr = dset.epr.SNR
   # Calculate the expected SNR based on the noise level and signal amplitude
    expected_snr = 20 * np.log10(1 / 0.04)  # Assuming signal amplitude is 1

    # Assert that the calculated SNR is within a tolerance of the expected SNR
    assert np.abs(dset_snr - expected_snr) / expected_snr < 0.1

def test_MeasurementTime(dataset_from_sequence):
    dataset = dataset_from_sequence
    dset = dataset.epr.correctphase
    assert dset.epr.MeasurementTime == 0.6

def test_merge_dataset():
    rng = np.random.default_rng(seed=0)
    def build_dataset(tmin):
        seq = T2RelaxationSequence(
            B=12200, LO=34.4, reptime=3e3, averages=1, shots=100, start=tmin)
        t, data = _simulate_T2(seq,0.2)
        data += rng.normal(0,0.04,data.shape) + 1j*rng.normal(0,0.04,data.shape)
        extra_params = {"nAvgs":1}
        dset = create_dataset_from_sequence(data,seq, extra_params)
        return dset
    
    dset1 = build_dataset(500)
    dset2 = build_dataset(dset1.tau[-1].values*1e3+10)

    dset = dset1.epr.merge(dset2)
    assert dset.data.shape[0] == dset1.data.shape[0] + dset2.data.shape[0]

