from pyepr import RectPulse, HSPulse, Detection
import numpy as np


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