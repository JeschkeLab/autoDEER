
Welcome to autoEPR & autoDEER!
====================================

.. caution:: 
    This documentation is a work in progress and not yet complete. 
    Please be patient and if anything is unclear, please contact the developers.

**autoEPR** is a free software package for the running of automated pulsed EPR 
experiments. It is designed to work with any arbitrary EPR sequence, including
shaped pulses, on any modern spectrometer. 

AutoEPR can be applied to many experiments in pulse EPR, currently, we are only
working on an implementation for DEER, autoDEER. 

**autoDEER** is a specific implementation of autoEPR for the purpose of push-button 
automatic Double Electron Electron Resonance (DEER) spectroscopy.

.. image:: source/images/autoDEER&autoEPR.svg

.. warning:: 
    autoEPR and autoDEER are actively developed software package, that is still very much a work in process. At this moment in time, we only recommend it is used by people who understand how
    the python code work and can debug when it goes wrong. Please consider this to be an alpha release, it is hoped that a more stable beta version can be released soon.
      
.. toctree::
    :maxdepth: 1
    :hidden:

    ./source/Install.rst
    ./source/autoDEER/index.rst
    ./source/autoEPR/index.rst
