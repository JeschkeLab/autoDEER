from autodeer.criteria import *
import time
import deerlab as dl
from autodeer.dataset import create_dataset_from_axes

def test_criteria():
    # Test initialization
    name = "test_criteria"
    description = "This is a test criteria"
    test = lambda x: x > 0
    criteria = Criteria(name, test, description)
    assert criteria.name == name
    assert criteria.description == description
    assert criteria.test == test

    # Test changing attributes
    new_name = "new_test_criteria"
    new_description = "This is a new test criteria"
    new_test = lambda x: x < 0
    criteria.name = new_name
    criteria.description = new_description
    criteria.test = new_test
    assert criteria.name == new_name
    assert criteria.description == new_description
    assert criteria.test == new_test

def test_TimeCriteria():
    now = time.time()
    criteria = TimeCriteria("Test Criteria", now, "")
    assert criteria.test(None) == True

def create_DEER_data(snr):
    t = np.linspace(0,3,50)
    r= np.linspace(2,6,20)
    Vmodel = dl.dipolarmodel(t,r,Pmodel=dl.dd_gauss)
    Vsim = Vmodel(mean=3,std=0.5,conc=150,scale=1,mod=0.7,reftime=0.4)

    params = {"tau1":400,"tau2":2600, "LO":34, "B":12220, "reptiem":3e3}
    Vsim += dl.whitegaussnoise(t,1/snr)
    return create_dataset_from_axes(Vsim, t, params, axes_labels=['t'])


def test_SNRCriteria():
    criteria = SNRCriteria(20)
    assert criteria.test(create_DEER_data(30)) == True
    assert criteria.test(create_DEER_data(15)) == False

def test_DEERCriteria():
    criteria = DEERCriteria("SPEED")
    assert criteria.test(create_DEER_data(50)) == True
    assert criteria.test(create_DEER_data(5)) == False