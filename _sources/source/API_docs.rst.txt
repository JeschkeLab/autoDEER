API doc
=================



Classes
----------------

Main Classes
~~~~~~~~~~~~

.. autoapisummary::
    
    autodeer.pulses.Pulse
    autodeer.criteria.Criteria
    autodeer.classes.Parameter
    autodeer.sequences.Sequence
    
Analysis Modules
~~~~~~~~~~~~~~~~

.. autoapisummary::

   autodeer.FieldSweep.FieldSweepAnalysis
   autodeer.ResPro.ResonatorProfileAnalysis
   autodeer.Relaxation.CarrPurcellAnalysis
   autodeer.Relaxation.ReptimeAnalysis
   autodeer.DEER_analysis.DEERanalysis

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

   autodeer.pulses.Pulse
   autodeer.pulses.Detection
   autodeer.pulses.Delay
   autodeer.pulses.RectPulse
   autodeer.pulses.GaussianPulse
   autodeer.pulses.HSPulse
   autodeer.pulses.ChirpPulse
   autodeer.pulses.SincPulse

Termination Criteria
~~~~~~~~~~~~~~~~~~~~

.. autoapisummary::

    autodeer.criteria.Criteria
    autodeer.criteria.TimeCriteria
    autodeer.criteria.SNRCriteria
    autodeer.criteria.DEERCriteria

Utilities
~~~~~~~~~

.. autoapisummary::

    autodeer.dataset.EPRAccessor
    autodeer.reporter.Reporter

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

    autodeer.DEER_analysis.DEERanalysis_plot
    autodeer.DEER_analysis.DEERanalysis_plot_pub
    autodeer.DEER_analysis.plot_overlap


Optimisation
~~~~~~~~~~~~

.. autoapisummary::

    autodeer.DEER_analysis.optimise_pulses
    autodeer.pulses.build_default_pulses
    autodeer.ResPro.optimise_spectra_position
I/O
~~~

.. autoapisummary::
    autodeer.tools.eprload
    autodeer.utils.save_file
    autodeer.dataset.create_dataset_from_sequence
    autodeer.dataset.create_dataset_from_axes
    autodeer.dataset.create_dataset_from_bruker

Utilities
~~~~~~~~~

.. autoapisummary::
    autodeer.utils.transpose_dict_of_list
    autodeer.utils.transpose_list_of_dicts
    autodeer.utils.round_step
    autodeer.DEER_analysis.normalise_01
    autodeer.utils.gcd
    autodeer.utils.sop
    autodeer.hardware.Bruker_tools.write_pulsespel_file