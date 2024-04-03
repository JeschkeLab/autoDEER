API doc
=================


Classes
----------------

Main Classes
~~~~~~~~~~~~

.. autoapisummary::
    
    autodeer.Pulse
    autodeer.Criteria
    autodeer.classes.Parameter
    autodeer.sequences.Sequence
    
Analysis Modules
~~~~~~~~~~~~~~~~

.. autoapisummary::

   autodeer.FieldSweepAnalysis
   autodeer.ResonatorProfileAnalysis
   autodeer.CarrPurcellAnalysis
   autodeer.ReptimeAnalysis
   autodeer.DEERanalysis

Sequences
~~~~~~~~~
.. _Sequences:
.. autoapisummary::

   autodeer.sequences.DEERSequence
   autodeer.sequences.HahnEchoSequence
   autodeer.sequences.T2RelaxationSequence
   autodeer.sequences.FieldSweepSequence
   autodeer.sequences.ReptimeScan
   autodeer.sequences.CarrPurcellSequence
   autodeer.sequences.RefocusedEcho2DSequence
   autodeer.sequences.ResonatorProfileSequence
   autodeer.sequences.TWTProfileSequence

Pulses
~~~~~~
.. _Pulses:

.. autoapisummary::

   autodeer.Pulse
   autodeer.Detection
   autodeer.Delay
   autodeer.RectPulse
   autodeer.GaussianPulse
   autodeer.HSPulse
   autodeer.ChirpPulse
   autodeer.SincPulse

Termination Criteria
~~~~~~~~~~~~~~~~~~~~

.. autoapisummary::

    autodeer.Criteria
    autodeer.TimeCriteria
    autodeer.SNRCriteria
    autodeer.DEERCriteria

Utilities
~~~~~~~~~

.. autoapisummary::

    autodeer.EPRAccessor
    autodeer.Reporter

Interfaces
~~~~~~~~~~
.. _Interfaces:

.. autoapisummary::
    autodeer.classes.Interface
    autodeer.hardware.Bruker_AWG.BrukerAWG
    autodeer.hardware.Bruker_MPFU.BrukerMPFU
    autodeer.hardware.XeprAPI_link.XeprAPILink
    autodeer.hardware.ETH_awg.ETH_awg_interface

Functions
----------------

Plotting
~~~~~~~~

.. autoapisummary::

    autodeer.DEERanalysis_plot
    autodeer.DEERanalysis_plot_pub
    autodeer.plot_overlap


Optimisation
~~~~~~~~~~~~

.. autoapisummary::

    autodeer.optimise_pulses
    autodeer.build_default_pulses
    autodeer.optimise_spectra_position
I/O
~~~

.. autoapisummary::
    autodeer.eprload
    autodeer.save_file
    autodeer.create_dataset_from_sequence
    autodeer.create_dataset_from_axes
    autodeer.create_dataset_from_bruker

Utilities
~~~~~~~~~

.. autoapisummary::
    autodeer.transpose_dict_of_list
    autodeer.transpose_list_of_dicts
    autodeer.round_step
    autodeer.normalise_01
    autodeer.gcd
    autodeer.sop
    autodeer.hardware.Bruker_tools.write_pulsespel_file