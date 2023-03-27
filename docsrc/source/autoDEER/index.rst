autoDEER
========================

AutoDEER has been designed to automate and simplify 
Double-Electron-Electron-Resonance (DEER), also known as 
Pulsed-Electron-Double-Resonance (PELDOR), experiments. DEER (or PELDOR) are 
a technique used to measure the distance between electron spin centres. This is 
normally applied to biological systems, by adding electron spin labels to specific
sites in a protein. 

Features of autoDEER include:

*   **Fully automated control:** A fully automated push-button mode controls the entire spectrometer, allowing the user to walk away. There also exists an advanced mode that allows the user to write their own automation scripts or to progress step by step. This gives an experienced user the ability to precisely control the measurement whilst still remaining easy and simple to use. 
*   **Support for a variety of DEER variants:** When run in push-button mode, the software will select a DEER sequence for the user. However, when using the advanced mode, the user can select a wide range of DEER variants, such as four-pulse, five-pulse, seven-pulse DEER and nDEER.
*   **Integrated analysis:** Both during and after a DEER measurement, the software actively processes and analyses the data using a tuned and optimised DeerLab-based implementation. [1]
*   **Parameter optimisation:** Sequence parameters are determined and optimised from the integrated analysis. This includes inter-pulse delays, pulse excitation bands, and pulse amplitudes.  
*   **Hardware independence:** autoDEER controls the spectrometer through a new Python-based EPR automation toolbox. This toolbox is spectrometer independent and supports both modern Bruker Elexsys-II (using the XeprAPI) and home-built spectrometers. Both traditional systems with multiple pulse forming units (MPFU) and modern AWG-based systems are supported. For Bruker systems, autoDEER can automatically write PulseSpel scripts, allowing the entire experiment to be written in Python.

.. toctree::
    :hidden:
    :maxdepth: 0
    :caption: Get started

    ./deer_intro.rst
    
.. toctree::
    :hidden:
    :maxdepth: 0
    :caption: Running an Experiment

    ./semi_autoDEER.rst
    ./fully_autoDEER.rst
    ./DEER_variants.rst

.. toctree::
    :hidden:
    :maxdepth: 0
    :caption: Post Processing

    ./analysis.rst

.. toctree::
    :hidden:
    :maxdepth: 0
    :caption: API Reference
    
    ./reference.rst


