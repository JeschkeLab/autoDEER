import os

# Set default values
waveform_precision = 2  # ns

def set_waveform_precision(precision):
    global waveform_precision
    waveform_precision = precision

def get_waveform_precision():
    global waveform_precision
    return waveform_precision
