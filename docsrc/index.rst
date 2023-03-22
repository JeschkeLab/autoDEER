
Welcome to autoDEER's documentation!
====================================

.. caution:: 
    This documentation is not yet complete. Please be patient and if anything is unclear, please contact the developers.

**autoEPR** is a free software package for the running of automated pulsed EPR 
experiments. It is designed to work with any arbitary EPR sequence, including
shaped pulses, on any modern spectrometer. 

**autoDEER** is a specific implementation of autoEPR for the purpose of push-button 
automatic Double Electron Electron Resonance (DEER) spectroscopy.

Features of autoDEER include:

*   **Fully automated control:** A fully automated push-button mode controls the entire spectrometer, allowing the user to walk away. There also exists an advanced mode that allows the user to write their own automation scripts or to progress step by step. This gives an experienced user the ability to precisely control the measurement whilst still remaining easy and simple to use. 
*	**Support for a variety of DEER variants:** When run in push-button mode, the software will select a DEER sequence for the user. However, when using advanced mode, the user can select a wide range of DEER variants, such as: four-pulse, five-pulse, seven-pulse DEER and nDEER.
*	**Integrated analysis:** Both during and after a DEER measurement, the software actively processes and analyses the data using a tuned and optimised DeerLab-based implementation. [1]
*	**Parameter optimisation:** Sequence parameters are determined and optimised from the integrated analysis. This includes: inter-pulse delays, pulse excitation bands, and pulse amplitudes.  
*	**Hardware independence:** autoDEER controls the spectrometer through a new Python-based EPR automation toolbox. This toolbox is spectrometer independent and supports both modern Bruker Elexsys-II (using the XeprAPI) and home-built spectrometers. Both traditional systems with multiple pulse forming units (MPFU) and modern AWG-based systems are supported. For Bruker systems, autoDEER can automatically write PulseSpel scripts, allowing the entire experiment to be written in Python.


.. warning:: 
    autoDeer is an actively developed software package, that is still very much a work in process. At this moment in time we only recomened it is used by people who understand how
    the python code work and can debug when it goes wrong. Please consider this to be an Alpha release, it is hoped that a more stable Beta version can be released soon.

      
.. toctree::
    :hidden:
    :maxdepth: 1

    ./source/Install.rst
    ./source/autoEPR/index.rst
    ./source/autoDEER/index.rst
    ./source/dev_guide/index.rst
