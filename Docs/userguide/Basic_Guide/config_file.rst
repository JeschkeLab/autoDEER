Configuration File
========================

Since this is a spectrometer independent software, the software must be
informed about the type of spectrometer being used and its specific limitations.
Even different generations of Bruker spectrometers can have very different hardware
limitations depending on its AWG.

The configuration file is writen as a YAML (Yet Another Markup Language) file,
this is an easy to write and easily computer and human readable.

An example: 

.. literalinclude:: ../../../test/test_data/test_Bruker_config.yaml
    :language: yaml

Options
------------
