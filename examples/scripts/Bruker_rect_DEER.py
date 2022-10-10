""" 
autoDEER MPFU
-------------

Running an automated rectangular DEER experiment on a BRUKER spectrometer using MPFU channels
---------------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------------

"""
import numpy as np
import autoDeer.hardware.xepr_experiments as exp
from autoDeer.ResPro import resonatorProfile
from autoDeer.FieldSweep import FieldSweep
from autoDeer.hardware import xepr_api as api
from autoDeer import std_deer_analysis
from scipy.io import savemat

gyro_e = 0.00281677
# Logs info
project = "autodeer"
sample = "MQ460_H:H_100uM_1_6mm"
folder = "2022_10_03_autoDEER_rect_FUS_attempt2"


# %%

xepr = api()  # This can only be run once per kernel
xepr.find_Xepr()
xepr.find_cur_exp()
xepr.find_hidden()

# This sets up the configuration file of spectrometer
if xepr.hidden['PlsPrgCalDbName'].value != 'Q_1101004':
    xepr.hidden['PlsPrgCalDbName'].value = 'Q_1101004'
    xepr.hidden['PlsPrgCalDbLoad']
    xepr.hidden['ApplyCfg']
    print("Changed config file")
else:
    print("Not changed config file")

# %%
# Measuring a first EDFS

# Set initial frequency guess
fc = 34.04
xepr.set_freq(fc)
xepr.set_field(fc/gyro_e)

# Tuning the necessary channels
tune = exp.MPFUtune(xepr, "Hahn", 16, 650)
tune.tune({'+<x>': "R+", '-<x>': "R-"}, tol=1, bounds=[20, 60])
# Change specjet base
d0 = (tune.calc_d0()//2)*2
print(d0)
# Run EDFS
exp.run_general(
    xepr,
    ["/PulseSpel/HUKA_DEER_AWG"],
    ["Field Sweep +<x>/-<x> pg=200ns", "Field sweep +<x>/-<x>"],
    {"PhaseCycle": False, "ReplaceMode": False},
    {"p0": 16, "p1": 16, "h": 25, "n": 1, "d0": 660}
    )

# Analyse EDFS

fs_data = xepr.acquire_scan()
fs = FieldSweep()
fs.import_from_dataclass(fs_data)
fs.find_max()

xepr.set_field(fs.max_field)
Bc = xepr.get_field()
fc = xepr.get_counterfreq()

fs.plot()
gyro_exp = fs.calc_gyro(xepr.get_counterfreq())
print(fs.gyro)
xepr.xepr_save(f"/home/xuser/xeprFiles/Data/HUKA/2022/{folder}/({project})_({sample})_(EDFS_Q)_init")
"""0.002816783144203381 """

# %%
# Perform resonator profile
t, nutation_data = exp.get_nutations(
    xepr,
    nu=[33.8, 34.4],
    field=[fc, Bc],
    step=20e-3)

rp = resonatorProfile()
res_prof = rp.calc_res_prof(nutation_data, dt=t[1]-t[0], freq_lim=[33.8, 34.4])
rp.autofit([33.8, 34.4])
rp.calculate_shape(res_prof, None)
rp.calc_IF_prof(2, rp.fc, type='fit')
rp.res_prof_plot([33.8, 34.4])

savedata = {"time": t, "Nutation Data": nutation_data}
savemat("/home/xuser/xeprFiles/Data/HUKA/2022/{folder}/({project})_({sample})_(ResPro)_init.mat",savedata)

# %%
# Tune ELDOR Channel
# DEER frequencies
fpump = fc + 0.01
fobs = fc - 0.04
xepr.set_ELDOR_freq(fpump)

# %%
# DEER time
# Setup DEER channels
xepr.set_freq(fobs)
xepr.set_field(fobs/gyro_exp)
tune = exp.MPFUtune(xepr, "Hahn", 32, 660)
tune.tune({'+<x>': "R+", '-<x>': "R-"}, tol=1, bounds=[20, 55])
tune = exp.MPFUtune(xepr, "Hahn", 16, 660)
tune.tune({'-<y>': "R+"}, tol=1, bounds=[30, 65])

xepr.set_ELDOR_freq(fpump)
xepr.set_field((fpump)/gyro_exp)

# %% 
# 4 pulse DEER
# 4p DEER
d0 = 660
tau1 = 800
tau2 = 2000
deadtime = tau1 - 80
scans = 500
step = 24
dim = np.floor((tau1 + tau2 - deadtime-100)/step)
exp.run_general(
    xepr,
    ["/PulseSpel/aD_DEER_MPFU"],
    ["4pDEER", "DEER run"],
    {"PhaseCycle": True, "ReplaceMode": False},
    {"p0": 32, "p1": 32, "p2": 32, "h": 20, "pg": 80, "n": scans, "d1": tau1,
     "d2": tau2, "d11": 200, "d3": deadtime, "d0": d0, "dim1": dim,
     "d30": step}, run=False)

xepr.run_exp()


# %% 
# Analyse DEER

run1 = xepr.acquire_scan()
fig, ROI, tau_max = std_deer_analysis(
    run1.axes / 1000, run1.data, tau1=tau1/1000, tau2=tau2/1000, 
    zerotime=0.08, precision="Normal", compactness=True)