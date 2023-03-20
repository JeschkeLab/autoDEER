# Here is a collection of pytests for the Xepr module. 
# This is very limited due to the obvious lack of a spectrometer.

import autoDeer.hardware.xepr_experiments as x_exp
import re
from autoDeer.hardware.XeprAPI_link import XeprApiLink 
import pytest

def test_change_dimensions():
    path = "autoDeer/PulseSpel/HUKA_DEER_AWG.exp"
    dim = 8
    new_dim = 256
    x_exp.change_dimensions(path, dim, new_dim)

    with open(path, 'r') as file:
        data = file.readlines()
    
    re_search = fr"dim{int(dim)}\s*s\[([0-9]+),*[0-9]*\]"

    for line in data:
        output = re.findall(re_search, line)
        if output != []:
            match = int(output[0])

    assert match == new_dim


def test_create_from_config():
    xepr = XeprApiLink("test/test_data/test_Bruker_config.yaml")

    assert xepr.AWG is False


def test_min_freq():
    xepr = XeprApiLink("test/test_data/test_Bruker_config.yaml")
    with pytest.raises(RuntimeError):
        xepr.set_freq(30)


def test_max_freq():
    xepr = XeprApiLink("test/test_data/test_Bruker_config.yaml")
    with pytest.raises(RuntimeError):
        xepr.set_freq(40)


def test_eldor_min_freq():
    xepr = XeprApiLink("test/test_data/test_Bruker_config.yaml")
    with pytest.raises(RuntimeError):
        xepr.set_ELDOR_freq(30)


def test_eldor_max_freq():
    xepr = XeprApiLink("test/test_data/test_Bruker_config.yaml")
    with pytest.raises(RuntimeError):
        xepr.set_ELDOR_freq(40)


def test_bridge_dt():
    xepr = XeprApiLink("test/test_data/test_Bruker_config.yaml")
    with pytest.raises(ValueError):
        xepr.set_PulseSpel_var("p1", 3)
