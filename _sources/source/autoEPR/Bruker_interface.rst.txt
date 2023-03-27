Bruker Interface
================

Since every Bruker Xepr spectrometer is slightly different and has a slightly 
different build specification a configuration file is required.

Nonetheless we do differentiate between an AWG-based spectrometer and a MPFU-based
spectrometer. These two generations of Bruker spectrometer represent a significant
change in both pulse and sequence definition. 

Configuration file
------------------

The configuration file is writen as a YAML (Yet Another Markup Language) file,
this is an easy to write whilst simultaneously both computer and human readable.

Here is a simple example: 

.. literalinclude:: ../../../test/test_data/test_Bruker_config.yaml
    :language: yaml


One important element of the configuration file is the frequency calibration.
This is only need if your Bruker spectrometer is using an analogue source and 
has not been upgraded to the newer digital source. When using analogue source,
you do not directly set the frequency but the potential difference (Voltage)
across the Gunn diode, this is done on a scale of 0-4095. For easy convertion a
polynomial fit of the Scale-Frequency graph is required. This must be at least
to the third order or preferably higher. In a DEER measurment frequencies have
to be very specific, down to the 1MHz. (i.e 0.03%).

Getting started
---------------
Before the Bruker interfaces can be imported the Bruker XeprAPI must be installed.
The recommended way of doing this is:

.. code-block:: bash

    python -m pip install "XeprAPI @git+https://github.com/OE-FET/XeprAPI"

Currently, the version on PyPI is incompatible with the latest numpy so we do not
recommended installing that version at this time. 


.. code-block:: python

    from autodeer.hardware import BrukerAWG # or BrukerMPFU

    interface = BrukerAWG("path/to/config/file")



