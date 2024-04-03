autoDEER GUI
============

1. Connecting to the spectrometer
---------------------------------
Before you can connect to a spectrometer you must first load a working directory.
It is in this working directory that all files will be saved.

File -> Load Folder

Next we need to load the configuration file. This is a ".yaml" file, 
that contains all the important spectrometer infomation, additionally this file 
can include preferences for autoDEER.

Spectrometer -> Load Config File

Finally, we must connect using:

Spectrometer -> Connect

If this is successful, the status-light should turn green. On Bruker spectrometer, 
this will take slightly longer as d0 needs to be estimated.


2. Starting an autoDEER Experiment
----------------------------------
All experiments require that the Maximum Measurment Time and Sample name are specified.  A fully automatic autoDEER experiment can be run by clicking the run autoDEER button. 

Advanced Mode
+++++++++++++
For advanced users there is the option to specify parameters. 
This could either be a specific sequence type or specify a inter-pulse delay.
Either a single inter-pulse delay is given and the others are found or all can be given. 

3. Viewing the analysis
----------------------------------

4. Printing a PDF lab-report
----------------------------------

Appendix
--------
