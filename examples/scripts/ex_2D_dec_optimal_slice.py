# %%
""" 
Optimal Slice Plots
-----------------------------------------------------------

Here we demonstrate how to get an optimal slice plot from some 2D Decoherence Data. This is useful for showing the improvements in signal gained by going from four
pulse DEER to 5 pulse DEER.

This data is from a doubly labelled WALP protein.

"""

import autoDeer as ad
from deerlab import deerload


t,V,params = deerload('../data/WALP23_H20_C7_C22_2Ddec_2Ddec_2hrs.DSC',full_output=True)

scans = int(params['DSL']['recorder']['NbScansDone'])
shots = int(params['DSL']['ftEpr']['ShotsPLoop'])
ShotRepTime = float(params['DSL']['ftEpr']['ShotRepTime'][:-3])

twoD_exp = ad.TwoD_Experiment()
twoD_exp.import_data(t,V,scans,shots,ShotRepTime)
fig = twoD_exp.optimal_slice_plot(norm='Max')
fig.show()

# %%
