Configuration File 
===================

As autoDEER is a spectrometer independent tool, it is necessary to specify the spectrometer configuration in a configuration file. 
The configuration file is a simple YAML file with the following format:

Overview
--------

* Spectrometer: 
  
  * Type: "Complete Spectrometer"
  
  * Manufacturer: "Bruker", "Bridge12" etc...
  
  * Model: "E580" etc...
  
  * Local Name: "My Spectrometer"

  * AWG: True/False

  * Bridge:
  
    * Min Freq: (in GHz)
  
    * Max Freq: (in GHz)
  
    * Sample Rate: (in GSa/s)
  
    * Det Freq: (in GSa/s)
  
    * Det Res: (in bits)
  
    * Power: (in W)

* Resonators:
  
  * (name of spectrometer):
    
    * Center Freq: (in GHz)
    
    * Q: (unitless)

* User Preferences:
  
  *  Initial Pulse Length: (in ns) [Used for initial experiments]
  

Examples
--------

BrukerMPFU
++++++++++
On Bruker spectrometers, the configuration files needs additional infomation.

* MPFU Channels: (A list of included MPFU Channels)

* Freq Cal: (A polynomical fit of the frequency calibration)
 
.. literalinclude:: ../../../config_files/BrukerMPFU.yaml
    :language: yaml
    :linenos:

ETHAWG
++++++++++

.. literalinclude:: ../../../config_files/ETHAWG.yaml
    :language: yaml
    :linenos:
