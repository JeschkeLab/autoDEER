# autoDeer
A Python package designed for the running of automated and optimised DEER/PELDOR experiments for pulsed EPR. This has been designed to be spectrometer independent and currently functions on both Bruker Elexsys-II and Andrin Doll style Spectrometers. 

## Features
- Automated tuning and pulse optimisation
- Automatic Distance characterisation and DEER parameter determination
- A generalised pythonic pulse EPR sequencer.
- Interfaces converting a general pulse sequence to either Matlab Structures or Bruker PulseSpel files 
- Easy development of new spectrometer interfaces
- 

## Dependencies
1) Numpy
2) Scipy
3) Matplotlib
4) pytest
5) h5py
8) DeerLab (https://jeschkelab.github.io/DeerLab/)

## Optinal Dependencies
1) Xepr API for python 3 (https://github.com/OE-FET/XeprAPI)
2) Matlab engine

## External Software Dependencies
1) Bruker Xepr
2) Matlab

## License

## Installation
Two main methods:
 - Either install from source by cloning this git repository
 - Install from a prebuilt release on the right hand side

Currently autoDEER is not published on any central repository such as pypi or conda. This will only happen after final publication.
