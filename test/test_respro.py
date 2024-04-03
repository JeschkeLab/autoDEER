from autodeer.ResPro import *
from autodeer.hardware.dummy import _similate_respro
from autodeer.sequences import ResonatorProfileSequence
from autodeer.dataset import create_dataset_from_sequence, create_dataset_from_axes
from autodeer.FieldSweep import create_Nmodel, FieldSweepAnalysis
import pytest

def test_ResonatorProfileAnalysis_from_sim():
    def lorenz_fcn(x, centre, sigma):
        y = (0.5*sigma)/((x-centre)**2 + (0.5*sigma)**2)
        return y
    
    mode = lambda x: lorenz_fcn(x, 34.0, 34.0/60)
    x = np.linspace(33,35)
    scale = 75/mode(x).max()
    mode_fun = lambda x: lorenz_fcn(x, 34.0, 34.0/60) * scale
    seq = ResonatorProfileSequence(B=12220,LO=34,reptime=3e3,averages=1,shots=50,fwidth=0.3)


    [tp_x, LO_axis], data = _similate_respro(seq,mode_fun)

    dset = create_dataset_from_sequence(data,seq)

    respro = ResonatorProfileAnalysis(dset
    )
    # respro.process_nutations(threshold=4)
    respro.fit()
    assert hasattr(respro,'results')
    assert hasattr(respro,'fc')
    assert hasattr(respro,'q')

    print(respro.fc)
    assert respro.fc == pytest.approx(34,abs=0.1)
    assert respro.q == pytest.approx(60,abs=2)

    fig = respro.plot()
    assert fig is not None



def test_ceil():
    assert ceil(3.5) == 4
    assert ceil(3.4) == 4
    assert ceil(3.0) == 3
    assert ceil(3.1) == 4
    assert ceil(3.9999999) == 4

    assert ceil(45,decimals=1) == 45
    assert ceil(45,decimals=-1) == 50

def test_floor():
    assert floor(3.5) == 3
    assert floor(3.4) == 3
    assert floor(3.0) == 3
    assert floor(3.1) == 3
    assert floor(3.9999999) == 3

    assert floor(45,decimals=1) == 45
    assert floor(45,decimals=-1) == 40

def test_calc_overlap():

    def gauss(x,centre,sigma):
        return np.exp(-((x-centre)/sigma)**2)
    
    gauss1 = lambda x: gauss(x,0,1)
    gauss2 = lambda x: gauss(x,0.5,1)

    x = np.linspace(-5,5,1000)
    overlap = calc_overlap(x,gauss1,gauss2)

    assert overlap == pytest.approx(0.352065)

def test_BSpline_extra():

    def gauss(x,centre,sigma):
        return np.exp(-((x-centre)/sigma)**2)
    
    x = np.linspace(-5,5,1000)
    y = gauss(x,0,1)

    spl = BSpline_extra((x,y,3))
    x_test = np.linspace(-10,10,100)
    y_test = spl(x_test)
    
    assert y_test[0] == pytest.approx(0,abs=1e-3)
    assert y_test[-1] == pytest.approx(0,abs=1e-3)
    assert y_test[50] == pytest.approx(1,abs=1e-1)

@pytest.fixture
def create_fake_respro():
    def lorenz_fcn(x, centre, sigma):
        y = (0.5*sigma)/((x-centre)**2 + (0.5*sigma)**2)
        return y
    
    mode = lambda x: lorenz_fcn(x, 34.0, 34.0/60)
    x = np.linspace(33,35)
    scale = 75/mode(x).max()
    mode_fun = lambda x: lorenz_fcn(x, 34.0, 34.0/60) * scale
    seq = ResonatorProfileSequence(B=12220,LO=34,reptime=3e3,averages=1,shots=50,fwidth=0.3)


    [tp_x, LO_axis], data = _similate_respro(seq,mode_fun)

    dset = create_dataset_from_sequence(data,seq)

    respro = ResonatorProfileAnalysis(
        nuts = dset.data.T,
        freqs = dset.LO.data,
        dt=2
    )
    respro.process_nutations(threshold=4)
    respro.fit()
    return respro

@pytest.fixture
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

def test_optimise_spectra_position(create_fake_fieldsweep,create_fake_respro):
    fsweep = create_fake_fieldsweep
    respro = create_fake_respro

    new_LO = optimise_spectra_position(respro,fsweep)

    assert new_LO == pytest.approx(-0.47,abs=0.1)  #TODO: check this value