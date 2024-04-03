The Algorithm
================

AutoDEER is based on a sophisticated algorithm to tune and optimize the DEER experiment.

The algorithm is based on the following steps:
1. A Echo Detected Field Sweep (EDFS) is performed to determine the precise spectrum of the sample, and for use in the optimization.
2. A nutation based resonator profile is measured to determine the resonators response function. This allows for pulses to be componsated for the resonator response.
3. Using the resonator profile and the EDFS optimal pulse shapes, bandwidths and lengths are determined.
4. A new EDFS is taken if the spectrometer frequency has been shifted.
5. A Carr-Purcell dynamical decoupling relaxation time measurement is performed to determine the relaxation times of the sample, and to calculate the optimal dipolar evolution time for the sample.
6. A short 2hr DEER experiment is performed to determine the region of interest (ROI) of the sample's distance distribution.
7. Using the ROI and the relaxation times, the optimal dipolar evolution time is calculated. And a longer time DEER experiment is performed.

This algorithm has been tested on a wide variety of samples and more detail can be seen here.


Tuning the pulses
-----------------


Echo Detected Field Sweep
-------------------------

Resonator Profile
-----------------

Relaxation measurements
-----------------------

