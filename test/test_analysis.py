from autodeer import eprload, FieldSweepAnalysis, CarrPurcellAnalysis, RectPulse, calc_optimal_deer_frqs, ChirpPulse
import numpy as np

def test_CarrPurcellAnalysis():
    dataset = eprload("test/test_data/test_CP.DSC")
    CP = CarrPurcellAnalysis(dataset)
    CP.fit("mono")
    CP.find_optimal(2*3600, 4, 40)
    pass


def test_FieldSweepAnalysis():
    dataset = eprload("test/test_data/test_FieldSweep.DSC")
    FS = FieldSweepAnalysis(dataset)
    FS.find_max()
    FS.calc_gyro(34.04)
    pass

def test_calc_optimal_deer_frqs_rect():
    dataset = eprload("test/test_data/test_FieldSweep.DSC")
    FS = FieldSweepAnalysis(dataset)
    FS.find_max()
    FS.calc_gyro(34.04)

    pump_pulse = RectPulse(tp=12,freq=0, flipangle=np.pi)
    exc_pulse = RectPulse(tp=12,freq=0, flipangle=np.pi/2)
    p_freq,e_freq = calc_optimal_deer_frqs(FS,pump_pulse,exc_pulse)

    assert np.abs((p_freq*1e3 - -15)/(p_freq*1e3)) < 0.15
    assert np.abs((e_freq*1e3 - -140)/(e_freq*1e3)) < 0.15

def test_calc_optimal_deer_frqs_rect_chirp():
    dataset = eprload("test/test_data/test_FieldSweep.DSC")
    FS = FieldSweepAnalysis(dataset)
    FS.find_max()
    FS.calc_gyro(34.04)

    pump_pulse = ChirpPulse(tp=12,init_freq=0,BW=0.15)
    exc_pulse = ChirpPulse(tp=12,init_freq=0, BW=0.05)
    p_freq,e_freq = calc_optimal_deer_frqs(FS,pump_pulse,exc_pulse)

    assert np.abs((p_freq*1e3 - -15)/(p_freq*1e3)) < 0.15
    assert np.abs((e_freq*1e3 - -140)/(e_freq*1e3)) < 0.15