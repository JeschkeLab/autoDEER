# autoDeer
A Python package designed for the running of automated and optimised DEER/PELDOR experiments for pulsed EPR. This has been designed to be spectrometer independent and currently functions on both Bruker Elexsys-II and Andrin Doll style Spectrometers. 

## Features
- Automated tuning and pulse optimisation
- Automatic Distance characterisation and DEER parameter determination
- A generalised pythonic pulse EPR sequencer.
- Interfaces converting a general pulse sequence to either Matlab Structures or Bruker PulseSpel files 

## Dependencies
1) Numpy
2) Scipy
3) Matplotlib
4) pytest
5) h5py
7) Xepr API for python 3 (https://github.com/OE-FET/XeprAPI)
8) DeerLab (https://jeschkelab.github.io/DeerLab/)

## External Software Dependencies
1) Bruker Xepr

## Troubleshooting
- <Can't find any open instances of Xepr>
    This has two main causes:
    1) Xepr API hasn't be started. Please do "Processing" -> "XeprAPI" -> "Enable Xepr API"
    2) The kernal has recently lost connection, and a background kernal is connected. Please open <htop> in terminal and
    kill the process of the other python kernal


