API Reference
====================


Core Functionality
+++++++++++++++++++++++
.. currentmodule:: autodeer.classes

.. autosummary::
    :toctree: _autosummary
        :maxdepth: 1
    :template: custom_class_template.rst
    :nosignatures:
    
    Dataset
    Interface
    Sequence
    Pulse
    Parameter
    Detection
    Delay



Default
++++++++++++++++++++
.. currentmodule:: autodeer.pulses

.. rubric:: Pulses
.. autosummary::
    :toctree: _autosummary
        :maxdepth: 1
    :template: custom_class_template.rst
    :nosignatures:
    
    RectPulse
    HSPulse
    ChirpPulse
    SincPulse



.. rubric:: Sequences
..  _standard_sequences: 
.. currentmodule:: autodeer.sequences

.. autosummary::
    :toctree: _autosummary
        :maxdepth: 1
    :template: custom_class_template.rst
    :nosignatures:
    
    HahnEchoSequence
    CarrPurcellSequence
    FieldSweepSequence
    DEERSequence

.. rubric:: Criteria


.. rubric:: Interfaces
.. currentmodule:: autodeer.hardware
    
.. autosummary::
    :toctree: _autosummary
        :maxdepth: 1
    :template: custom_class_template.rst
    :nosignatures:
    
    ETH_awg.ETH_awg_interface
    Bruker_MPFU.BrukerMPFU
    Bruker_AWG.BrukerAWG
    
    
    
