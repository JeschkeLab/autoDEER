autoEPR
========================

AutoEPR is a Python-based spectrometer-independent automation toolbox. It was originally designed as the back end for automating DEER experiments (autoDEER) but it is also readily expandable to a large variety of other pulse EPR experiments. Currently, it is still fully integrated into autoDEER this will likely change at a later date.

Features included in autoEPR:

* Automated Control: Fully automated experiments from pulse tunning and set-up experiments (Field-Sweeps, Resonator profiles, etc…) to final measurement and analysis. This can all be done such that the user only needs to press a single button and walk away. 
* Integrated Analysis: Whilst the experiment is running, the toolbox can actively process and analyse the data. This analysis can be used to either set future experimental parameters or to intelligently end the experiment once specific criteria have been satisfied. It is possible to expand to all Python-based data analysis packages.
* Generalised Sequences: Experimental sequences are written and designed in an easy-to-use object-based approach. Full support exists for both shaped and chirped pulses, including resonator and TWT compensation. Experiments and procedures written on one spectrometer are readily transferable to other spectrometers, even if they are of a different design or model.
* Hardware independence: This toolbox is spectrometer independent and supports both home-built spectrometers and modern Bruker Elexsys-II (using the XeprAPI). Both traditional systems with multiple pulse forming units (MPFU) and modern AWG-based systems are supported. For Bruker systems, autoEPR works by writing PulseSpel files itself. 


.. toctree::
    :hidden:
    :maxdepth: 1
    :caption: API Reference
    
    ./Getting_Started
    ./Sequencer
    ./Interface
    ./Reference