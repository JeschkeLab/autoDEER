from autodeer.tools import eprload
import numpy as np
from autodeer.utils import gcd

def test_eprload_mat():
    dataset = eprload("test/test_data/Matlab_file_test.mat")

    assert dataset.params["name"] == 'CHORUS'


def test_eprload_bruker():
    dataset = eprload("test/test_data/test_FieldSweep.DSC")
    assert dataset.B.min().data == 11993.65
    assert dataset.LO == 34.0645
    assert dataset.reptime == 5100.0