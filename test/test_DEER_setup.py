from autodeer.classes import *
from autodeer.dataset import *
from autodeer.sequences import *
from autodeer.pulses import *
from autodeer.DEER_analysis import *
import numpy as np
import pytest

from autodeer.hardware.dummy import _simulate_CP, _simulate_2D_relax, _simulate_T2, add_noise
from autodeer.Relaxation import CarrPurcellAnalysis, RefocusedEcho2DAnalysis, HahnEchoRelaxationAnalysis


def get_CPAnalysis(T_CP,V0=1.0):
    seq = CarrPurcellSequence(
        B=12220, LO=34.0, reptime=3e3,averages=1, shots=50,n=2, start=300,step=50, dim=200)
    x,V = _simulate_CP(seq,T_CP=T_CP)
    V *= V0
    V = add_noise(V,0.01)
    dataset = create_dataset_from_sequence(V,seq)
    dataset.attrs['nAvgs'] = 1
    CP = CarrPurcellAnalysis(dataset)
    CP.fit('auto')
    return CP

def get_TmAnalysis(Tm,V0=1.0):
    seq = T2RelaxationSequence(
        B=12220, LO=34.0, reptime=3e3,averages=1, shots=50, start=300,step=50, dim=200)
    x,V = _simulate_T2(seq,ESEEM_depth=0.1,Tm=Tm)
    V *= V0
    V = add_noise(V,0.01)
    dataset = create_dataset_from_sequence(V,seq)
    dataset.attrs['nAvgs'] = 1
    CP = HahnEchoRelaxationAnalysis(dataset)
    CP.fit('auto')
    return CP

def get_Ref2DAnalysis(sigma=0.8):
    seq = RefocusedEcho2DSequence(
        B=12220, LO=34.0, reptime=3e3,averages=1, shots=50, tau=10,)
    x,V = _simulate_2D_relax(seq,sigma=0.8)
    V = add_noise(V,0.01)
    dataset = create_dataset_from_sequence(V,seq)
    dataset.attrs['nAvgs'] = 1
    Ref2D = RefocusedEcho2DAnalysis(dataset)
    return Ref2D


def test_calc_deer_settings_Tm_CP_5pDEER():
    CP_dataset = get_CPAnalysis(4e3,V0=0.9)
    Tm_dataset = get_TmAnalysis(2e3)

    settings = calc_DEER_settings(CP_dataset,Tm_dataset,target_time=10,target_MNR=50)

    assert settings['ExpType'] == '5pDEER'
    assert settings['tau1'] == pytest.approx(6.4,abs=1e-2)
    assert settings['tau2'] == pytest.approx(6.4,abs=1e-2)
    assert settings['tau3'] == pytest.approx(0.3,abs=0.05)

def test_calc_deer_settings_Tm_CP_4pDEER():
    CP_dataset = get_CPAnalysis(0.9e3,V0=0.6)
    Tm_dataset = get_TmAnalysis(0.8e3)

    settings = calc_DEER_settings(CP_dataset,Tm_dataset,target_time=10,target_MNR=50)

    assert settings['ExpType'] == '4pDEER'
    assert settings['tau1'] == pytest.approx(0.4,abs=1e-2)
    assert settings['tau2'] == pytest.approx(1.9,abs=1e-2)

def test_calc_deer_settings_Ref2D_CP_5pDEER():
    CP_dataset = get_CPAnalysis(4e3,V0=0.9)
    Ref2D_dataset = get_Ref2DAnalysis()

    settings = calc_DEER_settings(CP_dataset,Ref2D_dataset,target_time=10,target_MNR=50)

    assert settings['ExpType'] == '5pDEER'
    assert settings['tau1'] == pytest.approx(6.4,abs=1e-2)
    assert settings['tau2'] == pytest.approx(6.4,abs=1e-2)
    assert settings['tau3'] == pytest.approx(0.3,abs=0.05)

def test_calc_deer_settings_Ref2D_CP_4pDEER():
    CP_dataset = get_CPAnalysis(0.9e3,V0=0.6)
    Ref2D_dataset = get_Ref2DAnalysis()

    settings = calc_DEER_settings(CP_dataset,Ref2D_dataset,target_time=10,target_MNR=50)

    assert settings['ExpType'] == '4pDEER'
    assert settings['tau1'] == pytest.approx(1.0,abs=1e-2)
    assert settings['tau2'] == pytest.approx(1.7,abs=1e-2)


