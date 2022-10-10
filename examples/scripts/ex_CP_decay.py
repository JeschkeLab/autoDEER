# %%
""" 
Analysing Carr-Purcel Decay Traces to extract optimal inter pulse delay
-----------------------------------------------------------

Here we look at how to quickly import a fine from bruker and calculate the 
optimal inter pulse delay.
"""

from autoDeer import Carr_Purcell


CP_example = Carr_Purcell()
CP_example.import_from_bruker("../data/CP_H20_example.DSC")

CP_example.fit()
plot = CP_example.plot()
optimal = CP_example.find_optimal(3000 * 1e-6, 20)
print(f"Optimal Dipolar Evolution Time: {optimal} us")
# %%
