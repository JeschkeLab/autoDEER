Configuration File
========================

Since this is a spectrometer independent software, the software must be
informed about the type of spectrometer being used and its specific limitations.
Even different generations of Bruker spectrometers can have very different hardware
limitations depending on its AWG.

The configuration file is writen as a YAML (Yet Another Markup Language) file,
this is an easy to write and easily computer and human readable.

An example: 

.. literalinclude:: ../../../test/test_data/test_Bruker_config.yaml
    :language: yaml

Options
------------

- Freq Cal:
    The frequency calibration is needed so that the correct frequencies can be set.
    Xepr does not allow you to set the bridge frequency directly and instead you
    set an integer of value 0-4095. There is a near linear relationship, between
    value and frequency however, a higher order polynomial fit is needed to get a 
    precise enough relationship. 

    The 'CalibrateFreq()' cmd can be used to generate the correct polynomial
    coefficients. This is best done when the spectrometer is still in operate mode,
    but with max attenuation and no pulses being used. The TWT can be off. 
