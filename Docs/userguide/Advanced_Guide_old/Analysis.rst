Analysis
===============

AutoDeer aims to analyse both decoherence and Deer experiemnts, and extract the necessary infomation to set further parameters. The analysis of Deer experiments
leans heavily on the work done by Fabregas-Iabenz, Stoll and Jeschke in DeerLab. [1]_ [2]_ For this reason DearLab v0.14 is required, and it is recomened that the user
reads and understands the DeerLab documentation. 

Decoherence experiments
-----------------------------

Decoherence experiemnts roughly fall in two catagories, 1 and 2 dimensional experiemnts. For 4p Deer, a 2D experiment is required as there is not yet a clear 
rule on finding optimal tau1 and tau2. [3]_ For 5 pulse deer and beyond, we are more interested in 1D experiemnts as we know how the delays relate to each other,
for this the work on Dynamical Decoupling is most applicable. [4]_


Deer experiemnts
----------------------

In autoDeer we have chosen to exclusively use DeerLab for the fitting and extraction of distance distributions for a few reasons. Primarily, because it is written 
in Python (like this package) but also as it can handle both 4 and 5 pulse experiments as well as maintaining a high degree of generailty. There is however a
consequence to this decision, and that comes at the cost of computational power. For this reason we recomened that  a moderatly powerfull modern computer is used,
this is often not the case with a Bruker spectrometer. To circumenvent this issue, we would like to add support for remote processing at a later date. This is not
possible currently and not an active priority, and subsequntly some features may not be possible on low powered computers.

AutoDeer also expands beyond DeerLab when it comes to parameter extraction, this package can find the region of intrest (ROI) and suggest appropriate time delays and
subsequent experiemnts


References
-------------------------
.. [1] Fábregas-Ibáñez, Luis, Stefan Stoll, and Gunnar Jeschke. “Compactness Regularisation in the Analysis of Dipolar EPR Spectroscopy.” Journal of Magnetic Resonance, 2022. 
        
.. [2] Fábregas Ibáñez, Luis, Gunnar Jeschke, and Stefan Stoll. “DeerLab: A Comprehensive Software Package for Analyzing Dipolar Electron Paramagnetic Resonance Spectroscopy Data.” Magnetic Resonance 1, no. 2 (October 1, 2020): 209–24. https://doi.org/10.5194/mr-1-209-2020.

.. [3] Bahrenberg, Thorsten, Samuel M. Jahn, Akiva Feintuch, Stefan Stoll, and Daniella Goldfarb. “The Decay of the Refocused Hahn Echo in Double Electron–Electron Resonance (DEER) Experiments.” Magnetic Resonance 2, no. 1 (April 16, 2021): 161–73. https://doi.org/10.5194/mr-2-161-2021.

.. [4] Soetbeer, Janne, Miriam Hülsmann, Adelheid Godt, Yevhen Polyhach, and Gunnar Jeschke. “Dynamical Decoupling of Nitroxides in O-Terphenyl: A Study of Temperature, Deuteration and Concentration Effects.” Physical Chemistry Chemical Physics 20, no. 3 (January 17, 2018): 1615–28. https://doi.org/10.1039/C7CP07074H.



.. .. toctree::
..     :maxdepth: 2
    
..     ./Decoherence_experiments