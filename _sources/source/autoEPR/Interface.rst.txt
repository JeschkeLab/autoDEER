Interfaces
==========

Currently, there are two main types of interface. 
1. A General Interface
2. A Bruker interface

A Bruker interface is still just a subclass of the general interface. However,
it has been specially designed to connect with Xepr. The general interface can
be expanded to support nearly all home-built spectrometers. 

.. image:: ../images/Spectrometer_interfaces.svg
    :align: center
    :width: 65%

Custom Interface
-----------------------
A key feature of this project is support for custom home-built spectrometers on
an equal footing to Bruker spectrometers. However, of course, we can not provide
native interfaces for every single spectrometer. Here is a short guide on how
to interface with a MATLAB-based spectrometer. This should give you a flavour 
of how to build your own.

All interfaces are a sub-class of Interface. This provides the minimum methods
required for an interface to work. 

Before starting to write an interface, you need to consider: How does my 
spectrometer read pulse sequences? and How do I convert the autoEPR pulse 
sequence to what is required?


Bruker Interface
-----------------------

.. caution:: 
    Due to stability issues development of the Bruker interface has paused,
    and is not recommended for use at this time. 
    Please contact the development team for further infomation. 

Most Bruker ELEXSYS-II spectrometers are fundamentally similar and are 
controlled by the Xepr. This means that we interface with all these
spectrometers in a similar fashion. Nonetheless, we still split them into 
two classes: MPFU-based and AWG-based. 

**Basic Approach**

Both of these interfaces are based upon the XeprAPI, which allows Python nearly
full control of Xepr with a couple of significant caveats. 

Most of the interface is done through the auto-generation of custom PulseSpel
scripts. This is done by the class PulseSpel and can be used separately if 
needed. 


**The need for calibration** 

Like all spectrometers, there is a need to perform a few calibration measurements
initially. On nearly-all spectrometers: a resonator profile and a TWT-profile 
are strongly recommended. However, on Bruker spectrometers, there are additional
requirements depending on your specific spectrometer. Many of these are caused
by non-digital components or non-linear components in the bridge and do not 
need to be done at the start of every measurement session.

+--------------------+---------------------------+  
| Calibration        | Required when?            |  
+====================+===========================+
| Resonator Profile  | Often                     | 
+--------------------+---------------------------+
| LO - Frequency     | Analogue Source           | 
+--------------------+---------------------------+

For more information please go to:

.. toctree::
    :maxdepth: 1
    
    ./Bruker_interface.rst