from autodeer.pulses import *
import pytest

def test_pulse_to_from_json():
    test_pulse = RectPulse(freq=0, tp=12, flipangle=np.pi/2)

    json_file = test_pulse._to_json()
    new_par = test_pulse._from_json(json_file)
    print(new_par)