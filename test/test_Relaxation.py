from autodeer.Relaxation import *
from autodeer import eprload
from autodeer.hardware.dummy import _simulate_CP, _simulate_reptimescan
from autodeer.sequences import DEERSequence, ReptimeScan
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
        tau1=3.0, tau2=3.0, tau3=0.3, dt=16, B=12220, LO=34.0, reptime=3e3, 
        averages=1, shots=50)

    seq.five_pulse(relaxation=True)
    x,V = _simulate_CP(seq)
    dataset = create_dataset_from_sequence(V,seq)
    CP = CarrPurcellAnalysis(dataset)
    CP.fit(fit_type)
    CP.find_optimal(2*3600, 4, 40, target_shrt=3e3, target_step=16)
    assert np.abs(CP.optimal - 5.9)/5.9 < 0.1

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