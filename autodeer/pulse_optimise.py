from pyepr import RectPulse, HSPulse, Detection
from pyepr.pulses import *
import numpy as np
from autodeer.DEER_analysis import optimise_pulses2, calc_est_modulation_depth, calc_functional


def build_default_pulses(AWG=True,SPFU=False,tp=12):
    """
    Create a dictionary of default pulses for a DEER experiment.

    Parameters
    ----------
    AWG : bool, optional
        If True, the pump pulse is a Hyperbolic Secant pulse, by default True
    SPFU : bool, optional
        If True, the pi pulse is double the length of the excitation pulse, by default False
    tp : int, optional
        Duration of the excitation pulse in nanoseconds, by default 12
    
    Returns
    -------
    dict
        Dictionary containing the excitation pulse, reference pulse, pump pulse and detection event.
    
    """

    if SPFU:
        tp_pi = tp*2
    else:
        tp_pi = tp

    exc_pulse = RectPulse(  
                    tp=tp, freq=0, flipangle=np.pi/2, scale=0
                )
    
    ref_pulse = RectPulse(  
                    tp=tp_pi, freq=0, flipangle=np.pi, scale=0
                        )
    
    if AWG:
        pump_pulse = HSPulse(  
                    tp=120, init_freq=-0.25, final_freq=-0.03, flipangle=np.pi, scale=0,
                        order1=6, order2=1, beta=10
                    )
        det_event = Detection(tp=512, freq=0)
    else:
        pump_pulse = RectPulse(
                    tp=tp_pi, freq=-0.07, flipangle=np.pi, scale=0)
        
        # det_event = Detection(tp=exc_pulse.tp * 2, freq=0)
        det_event = Detection(tp=128, freq=0)
        
    pulses = {'exc_pulse':exc_pulse, 'ref_pulse':ref_pulse, 'pump_pulse':pump_pulse, 'det_event':det_event}

    return pulses

def check_pulses_max_length(pulses, r_min=3.5, D = 52.04):
    """
    Checks all the pulses in either a dictionary or list to see if the maximum pulse length is exceeded.

    Equation:

    .. math::
        t_p = \\frac{1000}{4D/r_{min}^3}


    Parameters
    ----------
    pulses : dict or list
        Dictionary or list of pulses to check.
    r_min : float, optional
        Minimum distance to detect, this is used to calculate the maximum pulse length, by default 3.5nm
    D : float, optional
        Dipolar coupling constant in MHz nm^3, by default 52.04

    Returns
    -------
    bool
        True if the maximum pulse length is not exceeded, False otherwise.
    """
    max_length = 1000/(4*D/r_min**3) # maximum pulse length in ns

    if isinstance(pulses,dict):
        pulses = pulses.values()

    for pulse in pulses:
        if pulse.tp.value > max_length:
            return False
    return True
    

def create_pulses_rect(resonatorProfile, r_min=3.5, max_bandwidth=0.1, same_power=False):
    """
    Creates the optimal single frequency rectangular pulses for DEER using the resonator profile as a guide for the necessary power. r_max is used to set a maximum length.

    Parameters
    ----------
    resonatorProfile : pyepr.ResonatorProfile
        Resonator profile to be used as a guide for the pulse.
    r_min: float, optional
        Minimum distance to detect, this is used to calculate the maximum pulse length, by default 3.5nm
    max_bandwidth: float, optional
        Maximum bandwidth of the pulse in GHz, by default 0.1
    same_power: bool, optional
        If True, the pi refocusing pulses have the same power as the pi/2 excitation pulse otherwise they have the same bandwidth, by default False
        
    """
    D = 52.04 # MHz nm^3
    max_length = 1000/(4*D/r_min**3) # maximum pulse length in ns
    max_nu = 1/(2*max_length) # maximum nutation frequency in MHz

    min_length = 1/max_bandwidth

    resonator_bandwidth  = resonatorProfile.fc / resonatorProfile.q

    lower_3dB = resonatorProfile.model_func(resonatorProfile.fc-resonator_bandwidth)
    upper_3dB = resonatorProfile.model_func(resonatorProfile.fc+resonator_bandwidth)

    pulse_length = np.min([max_length,lower_3dB, upper_3dB])

    if pulse_length < min_length:
        pulse_length = min_length

    if same_power:
        pulse_length_pi2 = pulse_length/2
    else:
        pulse_length_pi2 = pulse_length
    
    exc_pulse = RectPulse(  
                    tp=pulse_length_pi2, freq=0, flipangle=np.pi/2, scale=0
                )
    
    ref_pulse = RectPulse(  
                    tp=pulse_length, freq=0, flipangle=np.pi, scale=0
                        )

    pump_pulse = RectPulse(
                    tp=pulse_length, freq=0, flipangle=np.pi, scale=0)
    
    det_event = Detection(tp=2*pulse_length, freq=0)

    pulses = {'exc_pulse':exc_pulse, 'ref_pulse':ref_pulse, 'pump_pulse':pump_pulse, 'det_event':det_event}

    return pulses 


def create_pulses_shape(
        resonatorProfile, spectrum, r_min=3.5, max_bandwidth=0.3, 
        ExcPulseShape=None, RefPulseShape=None, test_pulse_shapes=None,
        n_pump_pulses=1,verbosity=0):
    """
    Creates the optimal chirped shaped pulses for DEER using the resonator profile as a guide for the necessary power. r_max is used to set a maximum length.

    Algorithm
    ++++++++++
    1. Calculate the maximum pulse length based on r_min.
    2. Calculate the maximum nutation frequency based on the resoantor profile.
    3. Set the excitation pulse length to less than 90% of the maximum nutation frequency, and less that 50% of the bandwidth (resonator and spectral).
    4. Set the pump pulse bandwidth to the remaining bandwidth, either limited by the spectral bandwidth or the resonator.
    5. For each test pulse, default to [RectPulse,HSPulse,ChirpPulse], optimise the pulses and calculate the modulation depth.
    6. Select the pulse setup with the highest functional value.


    Parameters
    ----------
    resonatorProfile : pyepr.ResonatorProfile
        Resonator profile to be used as a guide for the pulse.
    spectrum: pyepr.FieldSweepAnalysis
    r_min: float, optional
        Minimum distance to detect, this is used to calculate the maximum pulse length, by default 3.5nm
    max_bandwidth: float, optional
        Maximum bandwidth of the pulse in GHz, by default 0.1
    test_pulse_shapes: list, optional
        List of pulse shapes to test, by default None. If None, the function will test RectPulse, HSPulse and ChirpPulse
    ExcPulseShape: pyepr.pulses.Pulse, optional
        Excitation pulse shape to use, by default None. If None, RectPulse is used
    RefPulseShape: pyepr.pulses.Pulse, optional
        Refocusing pulse shape to use, by default None. If None, RectPulse is used
    """

    D = 52.04 # MHz nm^3
    max_length = 1000/(4*D/r_min**3) # maximum pulse length in ns
    max_nu = 1/(2*max_length) # maximum nutation frequency in MHz

    if max_length > 256:
        max_length = 256 # Hard limit at 256ns

    if ExcPulseShape is None:
        ExcPulseShape = RectPulse
    if RefPulseShape is None:
        RefPulseShape = RectPulse
    if test_pulse_shapes is None:
        test_pulse_shapes = [RectPulse,HSPulse,ChirpPulse]

    _,_,spectum_bandwidth = spectrum.calculate_bandwidth(30) # 30dB bandwidth from spectrum
    resonator_bandwidth  = resonatorProfile.fc / resonatorProfile.q

    # Calculate the appropriate exc_pulse length
    max_nu = resonatorProfile.model.max()
    exc_pulse_length = 1/(2*max_nu*0.90)
    # Adjust the pulse length to be < 50% of availiable bandwidth
    bw_min_length = 1/(2*0.4*np.min([resonator_bandwidth,spectum_bandwidth]))
    if exc_pulse_length < bw_min_length:
        exc_pulse_length = bw_min_length
    
    if exc_pulse_length > max_length:
        exc_pulse_length = max_length

    
    
    exc_bandwidth = 1/exc_pulse_length # GHz
    pump_bandwidth = np.min([resonator_bandwidth,spectum_bandwidth]) - exc_bandwidth*0.5
    pump_bandwidth = np.around(pump_bandwidth,2) # round to 10 MHz

    if pump_bandwidth > max_bandwidth:
        pump_bandwidth = max_bandwidth
    
    PulseResults = {}
    FunctionalResults = {}
    ModDepthResults = {}
    for pump_pulse_type in test_pulse_shapes:
        exc_pulse = ExcPulseShape(tp=exc_pulse_length,flipangle=np.pi/2,scale=None)
        ref_pulse = RefPulseShape(tp=exc_pulse_length,flipangle=np.pi,scale=None)

        if pump_pulse_type == RectPulse:
            pump_pulse = pump_pulse_type(tp=exc_pulse_length,flipangle=np.pi,freq = -exc_bandwidth/2,scale=None)
        else:
            pump_pulse = pump_pulse_type(tp=max_length,flipangle=np.pi,init_freq=(-exc_bandwidth/2-pump_bandwidth),final_freq=-exc_bandwidth/2,scale=None)
            
        if verbosity > 0:
            print(f"Optimising pulses for {pump_pulse_type.__name__} with exc pulse length {exc_pulse_length:.2f} ns and pump bandwidth {pump_bandwidth:.3f} GHz")
        optimised_pulses,_ = optimise_pulses2(spectrum,pump_pulse,exc_pulse,ref_pulse,resonator=resonatorProfile,method = 'spectrum_shift',verbosity=verbosity,n_pump_pulses=n_pump_pulses)
        mod = calc_est_modulation_depth(spectrum,**optimised_pulses,respro=resonatorProfile,n_pump_pulses=n_pump_pulses)
        F = calc_functional(spectrum,**optimised_pulses,resonator=resonatorProfile,n_pump_pulses=n_pump_pulses)
        if verbosity > 0:
            print(f"Functional value for {pump_pulse_type.__name__}: {F:.4f}, Modulation depth: {mod:.4f}")
        
        PulseResults[pump_pulse_type] = optimised_pulses
        FunctionalResults[pump_pulse_type] = F
        ModDepthResults[pump_pulse_type] = mod

    best_pulse_type = max(FunctionalResults,key=FunctionalResults.get)

    results = PulseResults[best_pulse_type]
    det_event = Detection(tp=exc_pulse_length*2, freq=0)
    results['det_event'] = det_event
    return results


def create_pulse_from_dict(pulse_dict):
    """
    Create pulses from a dictionary of pulse parameters.

    Parameters
    ----------
    pulse_dict : dict
        Dictionary containing the pulse parameters.

    Returns
    -------
    dict
        Dictionary containing the created pulses.
    """
    pulses = {}
    for key, params in pulse_dict.items():
        pulse_type = params.pop('type')
        pulse_class = globals()[pulse_type]
        pulses[key] = pulse_class(**params)
    return pulses