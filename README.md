# autoDeer
A Python package designed for the running of automated and optimised DEER/PELDOR experiments for pulsed EPR, with minimal user interaction and training. This has been designed to be spectrometer independent and currently functions on both Bruker Elexsys-II and Andrin Doll style Spectrometers. 

## Features
- Automated tuning and pulse optimisation.
- Support for both rectangular and shaped pulses
- Easily runs: 4pulse, 5pulse or nDEER with equal ease.
- Automatic Distance characterisation and DEER parameter determination
- A generalised pythonic pulse EPR sequencer.
- Interfaces converting a general pulse sequence to either Matlab Structures or Bruker PulseSpel files 
- Easy development of new spectrometer interfaces
- 

## Dependencies
1) Numpy
2) Scipy
3) Matplotlib
4) [DeerLab](https://jeschkelab.github.io/DeerLab/) 
5) pyyaml
6) [XeprAPI](https://github.com/OE-FET/XeprAPI) [Bruker]
7) matlab-engine [Matlab]

## External Software Dependencies
1) Bruker Xepr [Bruker]
2) Matlab [Matlab based systems]

## License

## Installation
Two main methods:
 - Either install from source by cloning this git repository
 - Install from a prebuilt release on the right hand side

