Interfaces
================

Currently, there are two main types of interface. 
1. A general interface
2. A Bruker interface

A Bruker interface is still just a sub class of the general interface. However,
it has been specially designed to connect with Xepr. The general interface can
be expanded to support nearly all home-built spectrometers. 


Bruker Interface
-----------------------

Most Bruker ELEXSYS-II spectrometer are fundamentally similar and are 
controlled by the Xepr. This means that we interface with all these
spectrometers in a similar fashion. Nonetheless, we still split them into 
two classes: MPFU-based and AWG-based. 

**Basic Approach**

Both of these interfaces are based upon the XeprAPI, this allows Python nearly
full control of Xepr with a couple of significant caveats. 

Most of the interface is done through the auto-generation of custom PulseSpel
scripts. This is done by the class PulseSpel, and can be used seperately if 
needed. 


**The need for calibration** 

Like all spectrometers there is a need to perform a few calibration measurements
initially. On nearly-all spectrometer: a resonator profile and a TWT-profile 
are strongly recommended. However, on Bruker spectrometers there are additional
requirements depending on your specifici spectrometer. Many of these are caused
by non-digital components or non-linear components in the bridge, and only do 
not need to be repeated at the start of every measuremnt session.

+--------------------+---------------------------+  
| Calibration        | Required when?            |  
+====================+===========================+  
| Resonator Profile  | Always                    | 
+--------------------+---------------------------+  
| TWT Profile        | Always                    |   
+--------------------+---------------------------+  
| LO - Frequency     | Analogue Source           | 
+--------------------+---------------------------+  
| MPFU - Amplitude   | MPFU Channels (optional)  |  
+--------------------+---------------------------+  



Custom Interface
-----------------------
A key feature of this project is support for custom home-built spectrometer on
an equal footing to Bruker spectrometers. However, of course we can not provide
native interfaces for every single spectrometer. Here is a short guide on how
to interface with a MATLAB based spectrometer. This should give you a flavour 
of how to build your own.

All interfaces are a sub-class of Interface. This provides the minimum methods
required for an interface to work. 

Before starting to write an interface, you need to consider: How does my 
spectrometer read pulse sequences? and How do I convert the autoEPR pulse 
sequence to what is required?