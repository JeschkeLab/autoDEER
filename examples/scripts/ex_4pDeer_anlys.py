# %%
""" 
Analysing 4pDEER experiments to determine quality and the ROI
-----------------------------------------------------------

Here we look at how to quickly import a fine from bruker and abstract the MNR and the ROI.
"""


from autoDeer import std_4p_deer_analysis
from deerlab import deerload
import numpy as np

t,V = deerload("../data/WALP23_H20_C7_C22_2Ddec_DEER_d1_1_7us_d2_2400_2hrs_long.DSC")

ztime = t[np.argmax(V)] #If zerotime is not known it can be found

std_4p_deer_analysis(t,V,tau1=1.7,tau2=2.4,zerotime = ztime)


# %%
