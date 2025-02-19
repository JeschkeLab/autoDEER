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

## Requirements
AutoDEER is generally compataible with Windows, Mac and Linux and requires Python 3.8, 3.9, or 3.10.
The specific hardware implememntation may add additional limitation for example,
when using XeprAPI a Linux OS is needed.

### Dependencies
1) Numpy
2) Scipy
3) Matplotlib
4) [DeerLab](https://jeschkelab.github.io/DeerLab/) 
5) pyyaml
6) [XeprAPI](https://github.com/OE-FET/XeprAPI) [Bruker]
7) matlab-engine [Matlab]

### External Software Dependencies
1) Bruker Xepr [Bruker]
2) Matlab [Matlab based systems]

## Setup
At this time autoDEER is only avaliable from source. A packaged release will come later. 

## License
AutoDEER is licensed under the GNU GPLv3 public license, and is released without
warrenty or liability. Comercial use is allowed, however it is advised to contact the authors for support.

Copyright © 2021-2025: Hugo Karas, Stefan Stoll, and Gunnar Jeschke