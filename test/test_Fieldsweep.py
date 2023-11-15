
from autodeer.FieldSweep import *
from autodeer.dataset import *
import pytest

def test_create_Nmodel():
    model = create_Nmodel(34*1e3)

    assert model is not None
    assert hasattr(model, "Boffset")
    assert hasattr(model, "gy")
    assert hasattr(model, "gz")
    assert hasattr(model, "GB")
    assert hasattr(model, "az")
    assert hasattr(model, "axy")
    assert hasattr(model, "scale")

def test_create_Nmodel_with_params():
    model = create_Nmodel(34*1e3)
    B_axis = np.linspace(12100, 12300, 10)
    params = {"B":B_axis/1e1,"az":3.66,"axy":0.488,"gy":2.01,"gz":2.0,"GB":0.45,"scale":1, "Boffset":0.7}
    V = model(**params)
    ref_values = np.array([1.68557056e-11, 4.50164637e-01, 8.92010015e-01, 7.09120791e-01,
       4.82187836e-01, 3.09046641e-01, 1.07465055e-01, 1.33308381e-02,
       3.41919276e-13, 2.24934385e-34])
    assert np.allclose(V,ref_values) # refvalues from 15th Nov 2023

def test_FieldSweepAnalysis_from_Bruker():
    dset = create_dataset_from_bruker("test/test_data/test_FieldSweep.DTA")
    fsweep = FieldSweepAnalysis(dset)
    fsweep.calc_gyro()
    fsweep.fit()
    assert fsweep.gyro == 2.0023193043619996
    assert fsweep.results.gz == 2.0000000000000004

def test_erot():

    # Test with separate angles
    alpha = np.pi/4
    beta = np.pi/3
    gamma = np.pi/6
    R = erot(alpha, beta, gamma)
    expected_R = np.array([[-0.0474, 0.6597, -0.7500],
                           [-0.7891, 0.4356, 0.4330],
                           [0.6124, 0.6124, 0.500]])
    assert np.allclose(R, expected_R,rtol=1e-4,atol=1e-4)

    # Test with angles in an array
    angles = np.array([alpha, beta, gamma])
    R = erot(angles)
    assert np.allclose(R, expected_R,rtol=1e-4,atol=1e-4)

    # Test with invalid input
    with pytest.raises(ValueError):
        erot()
    with pytest.raises(ValueError):
        erot(alpha, beta)
    with pytest.raises(ValueError):
        erot(alpha, beta, gamma, "invalid_option")
    with pytest.raises(ValueError):
        erot(np.nan, beta, gamma)
    with pytest.raises(ValueError):
        erot(alpha, np.nan, gamma)
    with pytest.raises(ValueError):
        erot(alpha, beta, np.nan)

    # Test with option "rows"
    R1, R2, R3 = erot(alpha, beta, gamma, "rows")
    assert np.allclose(R1, expected_R[0, :],rtol=1e-4,atol=1e-4)
    assert np.allclose(R2, expected_R[1, :],rtol=1e-4,atol=1e-4)
    assert np.allclose(R3, expected_R[2, :],rtol=1e-4,atol=1e-4)

    # Test with option "cols"
    C1, C2, C3 = erot(alpha, beta, gamma, "cols")
    assert np.allclose(C1, expected_R[:, 0],rtol=1e-4,atol=1e-4)
    assert np.allclose(C2, expected_R[:, 1],rtol=1e-4,atol=1e-4)
    assert np.allclose(C3, expected_R[:, 2],rtol=1e-4,atol=1e-4)  

def test_eyekron():
    test_M = np.array([[1, 2], [3, 4]])
    new_M = eyekron(test_M)
    expected_M = np.array([[1, 2, 0, 0], [3, 4, 0, 0], [0, 0, 1, 2], [0, 0, 3, 4]])
    assert np.allclose(new_M, expected_M)

def test_kroneye():
    test_M = np.array([[1, 2], [3, 4]])
    new_M = kroneye(test_M)
    expected_M = np.array([[1, 0, 2, 0], [0, 1, 0, 2], [3, 0, 4, 0], [0, 3, 0, 4]])
    assert np.allclose(new_M, expected_M)

@pytest.fixture
def nitroxide_spinsystem():
    sys = SpinSystem([1/2],[1],[-0.0025 * 3.66 + 2.0175, 2.006, 2.003], [0.488*28.0328,0.488*28.0328,3.66*28.0328])
    sys.gn = np.array([0.4038])
    return sys

#TODO: Compare these tests direct withe easyspin
def test_ham(nitroxide_spinsystem):

    H = ham(nitroxide_spinsystem)
    assert H.shape == (6,6)

    assert np.allclose(H.diagonal(),np.array([51.300024+0j,0,-51.300024+0j,-51.300024+0j,0,51.300024+0j]))
    
def test_ham_eZ(nitroxide_spinsystem):
    [muxe,muye,muze] = ham_ez(nitroxide_spinsystem)
    assert np.allclose(muxe.diagonal(3), np.array([-14.05467926, -14.05467926, -14.05467926]))
    assert np.allclose(muxe.diagonal(-3), np.array([-14.05467926, -14.05467926, -14.05467926]))
    assert np.allclose(muye.diagonal(3), np.array([0+14.03823367*1j, 0+14.03823367*1j, 14.03823367*1j]))
    assert np.allclose(muye.diagonal(-3), np.array([-14.03823367*1j, -14.03823367*1j, -14.03823367*1j]))
    assert np.allclose(muze.diagonal(), np.array([-14.0172393+0.j, -14.0172393+0.j, -14.0172393+0.j,  14.0172393+0.j,
        14.0172393+0.j,  14.0172393+0.j]))
    
def test_ham_nz(nitroxide_spinsystem):
    [muxe,muye,muze] = ham_nz(nitroxide_spinsystem)
    assert np.allclose(muxe.diagonal(1), np.array([0.00217648+0j, 0.00217648+0j, 0.+0j, 0.00217648+0j,0.00217648+0j]))
    assert np.allclose(muxe.diagonal(-1), np.array([0.00217648+0.j, 0.00217648+0.j, 0.+0.j, 0.00217648+0.j,0.00217648+0.j]))
    assert np.allclose(muye.diagonal(1), np.array([0.-0.00217648j, 0.-0.00217648j, 0.+0.j, 0.-0.00217648j,0.-0.00217648j]))
    assert np.allclose(muye.diagonal(-1), np.array([0.+0.00217648j, 0.+0.00217648j, 0.+0.j, 0.+0.00217648j,0.+0.00217648j]))
    assert np.allclose(muze.diagonal(), np.array([ 0.003078+0.j,  0.+0.j, -0.003078+0.j,  0.003078+0.j,0.      +0.j, -0.003078+0.j]))
    
def test_resfields(nitroxide_spinsystem):
    Orientations = np.array([[0,0,1],[0,0,1],[0,0,1]])
    EigenFields, Intensities = resfields(nitroxide_spinsystem,Orientations,34.0)
    assert len(EigenFields) == 3
    assert len(Intensities) == 3

def test_build_spectrum(nitroxide_spinsystem):
    Bmin = 1200
    Bmax = 1250
    x,y = build_spectrum(nitroxide_spinsystem,34.0*1e3,(Bmin,Bmax))
    assert len(x) == len(y)
    assert x[0] == Bmin
    assert x[-1] == Bmax
    assert np.any(np.isnan(y)) == False


