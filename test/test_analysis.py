from autoDeer.Relaxation import Carr_Purcell
from autoDeer.FieldSweep import FieldSweep
from autoDeer import eprload


def test_CP():
    dataset = eprload("test/test_data/test_CP.DSC")
    CP = Carr_Purcell(dataset)
    CP.fit("mono")
    CP.find_optimal(2*3600, 4, 40)
    pass


def test_FieldSweep():
    dataset = eprload("test/test_data/test_FieldSweep.DSC")
    FS = FieldSweep(dataset)
    FS.find_max()
    FS.calc_gyro(34.04)
    pass
