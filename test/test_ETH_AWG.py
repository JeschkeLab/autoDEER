from autodeer.hardware.ETH_awg import ETH_awg_interface
from autodeer import HahnEchoSequence

# This tests without needing to be connected to the hardware 

def test_create_interface():

    interface = ETH_awg_interface(awg_freq=1.5, dig_rate=2, test_mode=True)
    assert interface.awg_freq == 1.5
    assert interface.dig_rate == 2
    assert interface.pulses == {}

def test_build_exp_struct():
    interface = ETH_awg_interface(awg_freq=1.5, dig_rate=2, test_mode=True)

    sequence = HahnEchoSequence(
        B=12200, LO=34.1, reptime=3e3, averages=1, shots=100)
    sequence.pulses[0].scale.value = 0.5
    sequence.pulses[1].scale.value = 0.5

    exp_struct = interface._build_exp_struct(sequence=sequence)

    assert exp_struct["LO"] == 34.1 - 1.5
    assert exp_struct["avgs"] == 1
    assert exp_struct["reptime"] == 3e3 * 1e3
    assert exp_struct["shots"] == 100
    assert exp_struct["B"] == 12200
    assert len(exp_struct["events"]) == 3

