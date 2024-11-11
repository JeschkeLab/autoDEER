Release Notes
=============

Version 0.9.0 (2024-11-11):
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


Version 0.8.0 (TBA):
+++++++++++++++++++++++++++

- Major Support Update for 2D Decoherence
  - `RefocusedEcho2DSequence` reforumlated to support match normal style
  - Added 2D Decoherence to the GUI
  - Added `RefocusedEcho2DAnalysis` class
- Improvements to Advanced User Mode
    - Added option to only measure up to RefocusedEcho2D
- New function `calc_deer_settings` to calculate DEER settings from the availiable relaxation data



Version 0.7.0 (2024-04-01):
+++++++++++++++++++++++++++

- Added Graphical User Interface (GUI)
- Major improvements to the automated algorithm and reliability
- Improvements to the PDF reports
- Updated Documentation and shift to autoapi based API docs
- Added initial support for Bruker AWG