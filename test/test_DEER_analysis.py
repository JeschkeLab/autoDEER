from autodeer.classes import *
from autodeer.dataset import *
from autodeer.sequences import *
from autodeer.pulses import *
from autodeer.DEER_analysis import *
from autodeer.FieldSweep import create_Nmodel, FieldSweepAnalysis
import numpy as np
from scipy.signal import hilbert
import pytest

def test_find_longest_pulse():
    seq = FieldSweepSequence(B=12200, Bwidth=500, LO=34.4, reptime=3e3, averages=1, shots=100)
    assert find_longest_pulse(seq) == 0.128


@pytest.mark.parametrize("compactness, model, ROI",[(True, None, False),
                                                    (False, None, False),
                                                    (False, dl.dd_gauss, False),
                                                    (False, dl.dd_gauss, True),
                                                    (False, dl.dd_gauss2, False),])
def testDEERanalysis_and_plot(compactness, model, ROI):
    
    
    dataset = create_dataset_from_bruker("test/test_data/test_DEER.DTA")
    dataset.coords['t'] = dataset.coords['t'] + 0.08 
    dataset.attrs['tau1'] = 2500
    dataset.attrs['tau2'] = 2500
    dataset.attrs['tau3'] = 200
    analysis = DEERanalysis(dataset,compactness=compactness, model=model, ROI=ROI)

    assert hasattr(analysis, "t")
    assert hasattr(analysis, "mask")
    assert hasattr(analysis, "r")
    assert hasattr(analysis, "Vexp")
    assert hasattr(analysis, "Vmodel")
    assert hasattr(analysis, "dataset")

    if ROI: 
        assert analysis.ROI is not None
        assert hasattr(analysis, "rec_tau_max")
        assert hasattr(analysis, "rec_dt")
    else:
        assert analysis.ROI is None
        assert not hasattr(analysis, "rec_tau_max")
        assert not hasattr(analysis, "rec_dt")
    
    if compactness:
        assert hasattr(analysis, "penweights")
    
    if model is None:
        assert hasattr(analysis, "P")
    elif model == dl.dd_gauss:
        assert hasattr(analysis, "mean")
    elif model == dl.dd_gauss2:
        assert hasattr(analysis, "mean1")
        assert hasattr(analysis, "mean2")

    assert analysis.stats["chi2red"] < 3 # Replace with data without crossing echoes

    # Test figure creation

    fig = DEERanalysis_plot(analysis, False, analysis.ROI)

    assert fig is not None


def test_IdentifyROI():
    r = np.linspace(1,6,100)
    P = dl.dd_gauss(r, mean=3, std=0.5)

    ROI = IdentifyROI(P, r)

    assert ROI[0] == 1.75
    assert ROI[1] == 4.25

def test_shift_pulse_freq_mono():
    pulse = RectPulse(tp=16,freq=0,flipangle=np.pi/2)
    pulse_shifted = shift_pulse_freq(pulse, 0.1)
    assert pulse_shifted.freq.value == 0.1

def test_shift_pulse_freq_chirp():
    pulse = ChirpPulse(tp=16, init_freq=0, BW=0.1, flipangle=np.pi/2)
    pulse_shifted = shift_pulse_freq(pulse, 0.1)
    assert pulse_shifted.init_freq.value == 0.1

    pulse = ChirpPulse(tp=16, final_freq=0, BW=0.1, flipangle=np.pi/2)
    pulse_shifted = shift_pulse_freq(pulse, 0.1)
    assert pulse_shifted.final_freq.value == 0.1

    pulse = ChirpPulse(tp=16, init_freq=0, final_freq=0.1, flipangle=np.pi/2)
    pulse_shifted = shift_pulse_freq(pulse, 0.1)
    assert pulse_shifted.init_freq.value == 0.1
    assert pulse_shifted.final_freq.value == 0.2

def test_normalise_01():
    X = np.linspace(-10,10,100)
    X_norm = normalise_01(X)
    assert X_norm.min() == 0
    assert X_norm.max() == 1

    X = np.linspace(2,10,100)
    X_norm = normalise_01(X)
    assert X_norm.min() == 0
    assert X_norm.max() == 1


def test_resample_and_shift_vector():
    # Test with a simple vector
    A = np.array([1, 2, 3, 4, 5])
    f = np.array([0, 1, 2, 3, 4])
    shift = 1
    expected_output = np.array([1, 1,2, 3, 4])
    assert np.allclose(resample_and_shift_vector(A, f, shift), expected_output)


def test_build_lowpass_butter_filter():
    filter_fun = build__lowpass_butter_filter(0.1)

    assert filter_fun is not None
    assert filter_fun(0) == 1
    assert np.isclose(filter_fun(0.1)[0], 0.157176)

def create_fake_fieldsweep():
    B_axis = np.linspace(12000, 12400, 100)
    Nmodel = create_Nmodel(34.0*1e3)
    params = {"B":B_axis/1e1,"az":3.66,"axy":0.488,"gy":2.01,"gz":2.0,"GB":0.45,"scale":1, "Boffset":0.7}
    V = Nmodel(**params)
    dset = create_dataset_from_axes(V, B_axis,params={"LO":34.0},axes_labels=['B'])
    fsweep = FieldSweepAnalysis(dset)
    fsweep.calc_gyro()
    fsweep.fit()
    return fsweep


def test_optimise_pulses():

    fsweep = create_fake_fieldsweep()
    pump_pulse = RectPulse(tp=16,freq=0,flipangle=np.pi)
    exc_pulse = RectPulse(tp=16,freq=0,flipangle=np.pi/2)
    ref_pulse = RectPulse(tp=16,freq=0,flipangle=np.pi)
    pump_pulse,exc_pulse,ref_pulse = optimise_pulses(fsweep,pump_pulse,exc_pulse,ref_pulse)

    assert np.abs(pump_pulse.freq.value - -0.042)/0.042 < 0.1
    assert np.abs(exc_pulse.freq.value - 0.012)/0.012 < 0.1
    assert np.abs(ref_pulse.freq.value - 0.012)/0.012 < 0.1

    fig = plot_overlap(fsweep, pump_pulse, exc_pulse, ref_pulse)
    assert fig is not None

from autodeer.hardware.dummy import _simulate_CP, _simulate_2D_relax
from autodeer.Relaxation import CarrPurcellAnalysis, RefocusedEcho2DAnalysis

def get_CPAnalysis():
    seq = CarrPurcellSequence(
        B=12220, LO=34.0, reptime=3e3,averages=1, shots=50, tau=10,n=2)
    x,V = _simulate_CP(seq)
    dataset = create_dataset_from_sequence(V,seq)
    dataset.attrs['nAvgs'] = 1
    CP = CarrPurcellAnalysis(dataset)
    CP.fit('double')
    return CP

def get_Ref2DAnalysis():
    seq = RefocusedEcho2DSequence(
        B=12220, LO=34.0, reptime=3e3,averages=1, shots=50, tau=10,)
    x,V = _simulate_2D_relax(seq)
    dataset = create_dataset_from_sequence(V,seq)
    dataset.attrs['nAvgs'] = 1
    Ref2D = RefocusedEcho2DAnalysis(dataset)
    return Ref2D


@pytest.mark.parametrize("exp",['auto','5pDEER','4pDEER'])
def test_calc_deer_settings(exp):
    CP_dataset = get_CPAnalysis()
    Ref2D_dataset = get_Ref2DAnalysis()

    settings = calc_deer_settings(exp,CP_dataset, Ref2D_dataset)

    if exp == '5pDEER' or exp == 'auto':
        assert settings['tau1'] == pytest.approx(2.478,rel=1e-3)
        assert settings['tau2'] == pytest.approx(2.478,rel=1e-3)
        assert settings['tau3'] == pytest.approx(0.3,abs=0.05)
        assert settings['AimTime'] == pytest.approx(2,rel=1e-3)
        assert settings['ExpType'] == '5pDEER'
    elif exp == '4pDEER':
        assert settings['tau1'] == pytest.approx(3.7,rel=1e-3)
        assert settings['tau2'] == pytest.approx(7.5,rel=1e-3)
        assert 'tau3' not in settings
        assert settings['AimTime'] == pytest.approx(2,rel=1e-3)
        assert settings['ExpType'] == '4pDEER'
     
    print(settings)