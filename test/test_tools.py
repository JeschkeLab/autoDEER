from autoDeer.tools import eprload
import numpy as np


def test_eprload_mat():
    dataset = eprload("test/test_data/Matlab_file_test.mat")

    assert dataset.params["name"] == 'CHORUS'


def test_eprload_bruker():
    dataset = eprload("test/test_data/test_DEER.DSC", full_output=True)

    assert type(dataset.axes) is np.ndarray
