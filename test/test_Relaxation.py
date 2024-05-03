from autodeer.Relaxation import *
from autodeer import eprload
from autodeer.hardware.dummy import _simulate_CP, _simulate_reptimescan, _simulate_2D_relax
from autodeer.sequences import DEERSequence, ReptimeScan, RefocusedEcho2DSequence
from autodeer.dataset import *
import pytest
from matplotlib.figure import Figure as mplFigure

# --------------------------- Test CarrPurcellAnalysis ------------------------

def test_CarrPurcellAnalysis_fromBruker():
    dataset = eprload("test/test_data/test_CP.DSC")
    CP = CarrPurcellAnalysis(dataset)
    CP.fit("mono")
    CP.find_optimal(2*3600, 4, 40, target_shrt=3e3, target_step=16)
    assert np.abs(CP.optimal - 3.2)/3.2 < 0.1

    fig = CP.plot()
    assert isinstance(fig,mplFigure)
    
@pytest.mark.parametrize("fit_type",["mono","double"])
def test_CarrPurcellAnalysis_dataset(fit_type):
    seq = DEERSequence(
        tau1=0.5, tau2=0.5, tau3=0.3, dt=16, B=12220, LO=34.0, reptime=3e3, 
        averages=1, shots=50)

    seq.five_pulse(relaxation=True, re_step=200)
    x,V = _simulate_CP(seq)
    dataset = create_dataset_from_sequence(V,seq)
    CP = CarrPurcellAnalysis(dataset)
    CP.fit(fit_type)
    CP.find_optimal(2*3600, 4, 40, target_shrt=3e3, target_step=16)
    assert np.abs(CP.optimal - 3.35)/3.35 < 0.1

    fig = CP.plot()
    assert isinstance(fig,mplFigure)


# --------------------------- Test ReptimeScanAnalysis ------------------------

def test_ReptimeAnalysis_from_dataset():
    seq = ReptimeScan(B=12220, LO=34.0, reptime_max=3e4, averages=1, shots=50)
    x,V = _simulate_reptimescan(seq)
    dataset = create_dataset_from_sequence(V,seq)
    RA = ReptimeAnalysis(dataset)
    RA.fit()
    optimal = RA.calc_optimal_reptime()

    assert np.abs(optimal - 3e3)/3e3 < 0.1
    fig = RA.plot()
    assert isinstance(fig,mplFigure)

# --------------------------- Test RefocusedEcho2DAnalysis ------------------------

def get_Ref2DAnalysis():
    seq = RefocusedEcho2DSequence(
        B=12220, LO=34.0, reptime=3e3,averages=1, shots=50, tau=10,)
    x,V = _simulate_2D_relax(seq)
    dataset = create_dataset_from_sequence(V,seq)
    dataset.attrs['nAvgs'] = 1
    Ref2D = RefocusedEcho2DAnalysis(dataset)
    return Ref2D

def test_RefocusedEcho2DAnalysis_from_dataset():
    Ref2D = get_Ref2DAnalysis()

    tau1, tau2 = Ref2D.find_optimal(type='4pDEER', SNR_target=20, target_time=2, target_step=15)
    assert tau1 == pytest.approx(4.0,rel=1e-3)
    assert tau2 == pytest.approx(7.9,rel=1e-3)

    tau1_2 = Ref2D.optimal_tau1(tau2)
    assert tau1_2 == pytest.approx(tau1,rel=1e-3)

    tau1, tau2 = Ref2D.find_optimal(type='5pDEER', SNR_target=20, target_time=2, target_step=15)
    assert tau1 == pytest.approx(6.9,rel=1e-3)
    assert tau2 == pytest.approx(6.9,rel=1e-3)


    

    
