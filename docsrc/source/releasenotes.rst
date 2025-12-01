Release Notes
=============
Version 1.0.1 (2025-11-23):
++++++++++++++++++++++++++++
- Fixed the definition of the exciation probabiklity profile for pi/2 pulses to be the transverse component (sqrt(Mx^2 + My^2)) rather than (-Mz + 1), which was incorrect.
- 

Version 1.0.0 (2025-09-12):
++++++++++++++++++++++++++++
- All references to `LO` have been changed to `freq` in the frequency object and related.
- Moved to a Material UI theme for GUI.
- Improved pulse optimisation and optimisation.
- Added 1D refocused experiments.
- Fully tested protocol, details in the paper.


Version 0.10.0 (2024-12-01):
+++++++++++++++++++++++++++

- Split of PyEPR into a separate package
  - autoDEER now uses PyEPR as a dependency
  - Many functions and classes have been moved to PyEPR
  - Some paths have been changed to reflect the new package structure
- All sequences now have a simulate method that is used by the dummy interface.
- Move to Poetry for package management, and the default installation method
- Version numbering now stored in pyproject.toml


Version 0.9.0 (2024-11-11) (with PyEPR):
+++++++++++++++++++++++++++

- Faster implmenetation of excite profile
- Implementation of DeerLab based fitting into relaxation traces
- New implementation of correction factor
- Updated MNR and SNR parameters
- Improved versioning (including git branch for dev versions)
- Added a button for keep running (not to automatically stop)
- Implementation of a global waveform precision
- Added autoextending relxation data traces
- New save folder for raw data folder
- Improved calc_DEER_delays function, with new plot


Version 0.8.0 (TBA) (with PyEPR):
+++++++++++++++++++++++++++

- Major Support Update for 2D Decoherence
  - `RefocusedEcho2DSequence` reforumlated to support match normal style
  - Added 2D Decoherence to the GUI
  - Added `RefocusedEcho2DAnalysis` class
- Improvements to Advanced User Mode
    - Added option to only measure up to RefocusedEcho2D
- New function `calc_deer_settings` to calculate DEER settings from the availiable relaxation data



Version 0.7.0 (2024-04-01) (with PyEPR):
+++++++++++++++++++++++++++

- Added Graphical User Interface (GUI)
- Major improvements to the automated algorithm and reliability
- Improvements to the PDF reports
- Updated Documentation and shift to autoapi based API docs
- Added initial support for Bruker AWG