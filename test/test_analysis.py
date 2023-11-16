from autodeer import eprload, FieldSweepAnalysis, CarrPurcellAnalysis, RectPulse, optimise_pulses, ChirpPulse
import numpy as np

def test_CarrPurcellAnalysis():
    dataset = eprload("test/test_data/test_CP.DSC")
    CP = CarrPurcellAnalysis(dataset)
    CP.fit("mono")
    CP.find_optimal(2*3600, 4, 40, target_shrt=3e3, target_step=16)
    pass


def test_FieldSweepAnalysis():
    dataset = eprload("test/test_data/test_FieldSweep.DSC")
    FS = FieldSweepAnalysis(dataset)
    FS.find_max()
    FS.calc_gyro(34.04)
    FS.fit()
    pass

def test_calc_optimal_deer_frqs_rect():
    dataset = eprload("test/test_data/test_FieldSweep.DSC")
    FS = FieldSweepAnalysis(dataset)
    FS.find_max()
    FS.calc_gyro(34.04)
    FS.fit()

    pump_pulse = RectPulse(tp=12,freq=0, flipangle=np.pi)
    exc_pulse = RectPulse(tp=12,freq=0, flipangle=np.pi/2)
    pump_pulse, exc_pulse, ref_pulse = optimise_pulses(FS, pump_pulse, exc_pulse, exc_pulse)
    p_freq = pump_pulse.freq.value *1e3
    e_freq = exc_pulse.freq.value *1e3
    assert np.abs(p_freq - -80) < 5
    assert np.abs(e_freq - -6) < 5

def test_calc_optimal_deer_frqs_rect_chirp():
    dataset = eprload("test/test_data/test_FieldSweep.DSC")
    FS = FieldSweepAnalysis(dataset)
    FS.find_max()
    FS.calc_gyro(34.04)
    FS.fit()

    pump_pulse = ChirpPulse(tp=12,init_freq=0,BW=0.15,flipangle=np.pi)
    exc_pulse = ChirpPulse(tp=12,init_freq=0, BW=0.05,flipangle=np.pi/2)
    pump_pulse, exc_pulse, ref_pulse = optimise_pulses(FS, pump_pulse, exc_pulse, exc_pulse)
    
    p_freq = pump_pulse.init_freq.value *1e3
    e_freq = exc_pulse.init_freq.value *1e3
    
    assert np.abs(p_freq - -160) < 5
    assert np.abs(e_freq - -30) < 5