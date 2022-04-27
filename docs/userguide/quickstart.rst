Quickstart Guide
==================

AutoDeer aim is to provide an as automated as possible experiance for Double Electron Electron Resonance (DEER). 
It has been written to be used with three types of spectrometers (Bruker, Bruker-Hybrid and Homebuilt).

Spectrometer types
-------------------

*Bruker*

This uses the Bruker Xepr Python API to control an E580 Spectrometer. The Spectrometer must have the digital upgrade that enables the API, and have 4 MPFU channels.
There is currently no support for Bruker AWG, as we do not have access to one. 

*Bruker-Hybrid*

This is designed to work with a Bruker E580 Spectrometer where a Keysight M8190A AWG is used as an incoherent source for probe pulses. 
The AWG is triggered by the ELDOR channel on the E580 Spectrometer.

*Homebuilt*

Currently there is no control for Homebuilt spectrometers, due to the wide range of implementations. 
The automated analysis does work when using an old-style ".mat" file to save data.


Features:
----------

*Control*

Currently, autoDeer can control Bruker E580 Spectrometers and Keysight M8190A AWG, within a limited set of scenarios.
All of this control is still very buggy and unreliable and is *not* yet recomened for regular use. This will hopefully 
change withing the next few months.

*Analysis*

The Second part of an automated spectrometer is analysing data and choosing parameters. This can take the form of a decohernce
experiment or of a 4/5 pulse Deer experiment. This is much more stable and can be used. 


Interfaces:
-----------

There are two main ways of interfacing with this code.

* Jupyter Notebooks (recomened)
* Scripts
* GUI (comming soon)

Using Jupyter Notebooks is recomened at the moment and example notebooks are provided. However, on some setups Bruker's API can be more prone t
segmentaion faults when using Jupyter Notebooks, in this case running it as a script on in console is recomened, albeit annoying. 

