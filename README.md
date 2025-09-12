# autoDEER

## About
autoDEER is a Python based software package for running automated DEER experiments on a variety of hardware developed by the Jeschke Lab @ ETH. It is built on top of PyEPR, a general purpose EPR data processing and simulation library also developed by the Jeschke Lab @ ETH.


## Features
- Automated DEER experiments with real-time data analysis and adaptive experiment control
- Support for multiple hardware platforms (currently Bruker ElexSys and home-built spectrometers)
- Modular and extensible design for easy integration of new hardware and experiment types
- User-friendly interface for setting up and monitoring experiments
- Comprehensive documentation and examples to get started quickly
- Advanced users can run autoDEER from scripts for full control and customization

## Installation
Due to the variety of hardware and software dependencies, it is recommended to install autoDEER from source. At the moment, no pre-compiled packages are available.
Please refer to the [Installation Guide](https://jeschkelab.github.io/autoDEER/source/Install.html) in the documentation for detailed instructions.

### Dependencies
Non hardware specific dependencies:
- numpy
- scipy
- matplotlib
- [PyEPR](https://jeschkelab.github.io/PyEPR/)
- [DeerLab](https://jeschkelab.github.io/DeerLab/)
- pyyaml
- xarray
- h5netcdf
- toml
- numba
- PyQt6
- qt-material
- threadpoolctl
- quadprog
- reportlab

## User Guide
Please refer to the [User Guide](https://jeschkelab.github.io/autoDEER/gui_guide.html) in the documentation for detailed instructions on how to use autoDEER.

## Citing autoDEER
The paper associated with autoDEER has been submitted. Once published, it will be linked here. It is kindly requested that this paper is cited when using autoDEER in your research. Additionally, please cite DeerLab when using autoDEER.

Citing academic software is important as it helps to ensure the long-term sustainability of the software, and allows the developers to track the impact of their work and secure future funding. It also helps to provide credit to the developers for their hard work.

## Contributing
Contributions to autoDEER are welcome! If you have discovered an issue or have a feature request, please open an issue on the GitHub repository. If you would like to contribute code, please fork the repository and submit a pull request. If you have any questions or need help, please open an issue or contact the authors.

## License
autoDEER is licensed under the GNU GPLv3 public license, and is released without
warranty or liability. Commercial use is allowed, however it is advised to contact the authors for support.

Copyright © 2021-2025: Hugo Karas, Stefan Stoll and Gunnar Jeschke




