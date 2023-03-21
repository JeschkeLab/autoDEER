API Reference
====================


Core Functionality
+++++++++++++++++++++++
.. currentmodule:: autodeer.openepr

.. autosummary::
    :toctree: _autosummary
        :maxdepth: 1
    :template: custom_class_template.rst
    :nosignatures:
    
    dataset
    Interface
    Sequence
    Pulse
    Parameter


Default
++++++++++++++++++++

.. rubric:: Pulses
.. autosummary::
    :toctree: _autosummary
        :maxdepth: 1
    :template: custom_class_template.rst
    :nosignatures:
    
    Detection
    Delay
    RectPulse
    HSPulse
    ChirpPulse
    SincPulse
    ChorusPulse



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
    
    ETH_awg_interface
    BrukerMPFU
    BrukerAWG
    