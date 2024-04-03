GUI Interface
=============

The GUI interface is the recommened way to use autoDEER.
It is a simple porgram that allows you to run and analyse autoDEER experiments,
with no programing knownledge required. Additionally, it still allows minor
customisation of the experiment, to suit a users specific preferences. 

The GUI also has advanced features such as the ability to generate a PDF report. 

Connecting to a spectrometer
----------------------------

1. Open the GUI interface by running:
2. Click the "Spectrometer -> Load Config File" button and choose your spectrometer config file.
3. Click the "Spectrometer -> Connect" button.
4. A green light should appear in the top right corner of the GUI, indicating that the spectrometer is connected.
5. A folder must be opened to save the data to. This can be done by clicking the "File -> Open Folder" button. The path will be displayed below the 'spectrometer connected' light.

Starting an experiment
----------------------

Fully automatic experiment (recommended)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. On the first tab screen
2. Enter the maximum measurement time in the "Max Time" box.
3. Enter your sample name in the "Sample Name" box.
4. Enter the measurement temperature in the "Temperature" box.
5. Select the resonator you wish to use in the "Resonator" box, and adjust the "Frequency" box if necessary.
This frequency is the initial guess for the resonator frequency, it does not need to be the center but must be within the range of the resonator. 

6. Click the "Run Fully Automatic DEER" button.

Advanced Mode
~~~~~~~~~~~~~
1. Complete all the fields in the "Fully Automatic DEER" section, except don't press the button.
2. Input your preferences on the right hand side of the GUI. Any fields left blank will still be calculated automatically. 
3. Click the "Run Advanced Mode" button.

During an Experiment
--------------------
The experiment will now run automatically. As the data is recieved and processed, it will be displayed in the GUI.
At the bottom of the GUI a text box will display the current status of the experiment.

After an Experiment
-------------------
Once the experiment has finished or during an experiment, a PDF report can be now be saved. "File -> Save Report" will save a PDF report, in a chosen folder.
The report will contain all the data and graphs from the experiment, as well as a summary of the experiment parameters. This is somewhat akin to a lab book entry.
It is recommended that this along with all saved data is kept for future reference.
Data will be saved in both the spectrometers native format and as a '.h5' file. The '.h5' file is a HDF5 file, which is a standardised data format and used widely in the scientific community. It can be easily read by python and matlab, and is recommended for use with autoDEER.
 