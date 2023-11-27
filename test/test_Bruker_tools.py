from autodeer.hardware.Bruker_tools import *
from autodeer.hardware.Bruker_MPFU import _MPFU_channels

from autodeer.sequences import *
import numpy as np
import pytest




@pytest.mark.parametrize("seq", ['deer','fieldsweep','respro', 'reptime'])
def test_write_pulsespel_file_MPFU(seq):
    if seq == 'deer':
        seq = DEERSequence(tau1=3.0, tau2=3.0, tau3=0.2, dt=12, B=12200, LO=34, reptime=3e3, averages=100, shots=100)
        seq.five_pulse()
        seq.select_pcyc("DC")
    elif seq == 'fieldsweep':
        seq = FieldSweepSequence(
            B=12200,Bwidth=500, LO=34.0, reptime=3e3, averages=1, shots=10)
    elif seq == 'respro':
        seq = ResonatorProfileSequence(
            B=12200, LO=34.0, reptime=3e3, averages=1, shots=10)
    elif seq == 'reptime':
        seq = ReptimeScan(B=12200, LO=34.0,reptime=3e3, reptime_max=3e4, averages=1, shots=10)
    
    channels = _MPFU_channels(seq)
    MPFU = ["+<x>","-<x>","+<y>","-<y>"]
    def_file, exp_file = write_pulsespel_file(seq,AWG=False,MPFU=MPFU)
    if seq == 'deer' or seq == 'respro':
        assert 'sweep' in exp_file
    elif seq == 'fieldsweep':
        assert 'bsweep' in exp_file
    elif seq == 'reptime':
        assert 'sweep' not in exp_file

    print(exp_file)

def test_build_unique_progtable():
    seq = ResonatorProfileSequence(
            B=12200, LO=34.0, reptime=3e3, averages=1, shots=10)
    uprogtable = build_unique_progtable(seq)
    print(uprogtable)
    n_dim = len(uprogtable)
    assert n_dim == 2
    # Check variable order is correct
    assert uprogtable[0]['variables'][0]['variable'] == (0,'tp')
    assert uprogtable[1]['variables'][0]['variable'] == (None,'B')
    assert uprogtable[1]['variables'][1]['variable'] == (None,'LO')
    # Check axis is correct to sequence
    assert uprogtable[0]['axis']['start'] ==  seq.pulses[0].tp.axis[0]['axis'].min()
    assert uprogtable[0]['axis']['step'] ==  seq.pulses[0].tp.axis[0]['axis'][1]- seq.pulses[0].tp.axis[0]['axis'][0]
    assert uprogtable[0]['axis']['dim'] ==  seq.pulses[0].tp.axis[0]['axis'].shape[0]


def test_PulseSpel_ResonatorProfile():

    pass
def test_PulseSpel_TWTProfile():
    pass