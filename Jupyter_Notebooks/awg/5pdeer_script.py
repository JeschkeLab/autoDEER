import numpy as np
import matplotlib.pyplot as plt
import autoDeer.hardware.keysight_awg as ksawg
import autoDeer.hardware.pulses as pulses
import autoDeer.hardware.xepr_experiments as exp
import autoDeer.hardware.awg_experiments as awg_exp
from autoDeer.ResPro import resonatorProfile
from autoDeer.FieldSweep import FieldSweep
from autoDeer.hardware import xepr_api as api
import autoDeer.Param_Optimization as po
import time

awg = ksawg.Interface()
awg.open_con('129.132.218.87')
awg._inst_init()
print(f"Checking Sequencing mode:{awg.getFunctionMode(3)}")
sampling_freq = 12e9
sampling_period = 1/sampling_freq
grad = 64

xepr=api() # This can only be run once per kernel
xepr.find_Xepr()
xepr.find_cur_exp()
xepr.find_hidden()

fc = 34.447
Bc = 12227.90
xepr.set_field(Bc)
xepr.set_freq(fc)

exp.run_general(xepr,
    ["/PulseSpel/HUKA_DEER_AWG"],
    ["Field Sweep +<x>/-<x> pg=200ns","Field sweep +<x>/-<x>"],
    {"PhaseCycle":False,"ReplaceMode":False},
    {"p0":16,"p1":16,"h":25,"n":1,"d0":600})

time.sleep(10)

fs_data = xepr.acquire_scan()
fs = FieldSweep()
fs.import_from_dataclass(fs_data)
fs.find_max()

xepr.set_field(fs.max_field)
Bc = xepr.get_field()
fc = xepr.get_counterfreq()
print(f"max field = {fs.max_field:.1f}")
# fs.plot()

awg_exp.sequence_nutation(awg,1,1,128,20)

awg.Abort()
awg.StartSigGen(3)

exp.run_general(xepr,
    ["/PulseSpel/HUKA_DEER_AWG"],
    ["AWG seq inv (p2)","AWG +-<x> obs"],
    {"PhaseCycle":True},
    {"p0":16,"p1":16,"h":10,"n":1,"d0":580},False)

time.sleep(2)
print('Starting exp')
xepr.run_exp()
time.sleep(2)
print('Acquiring Scan')
run1 = xepr.acquire_scan()

print('Starting exp')
xepr.run_exp()
time.sleep(2)
print('Acquiring Scan')
run2 = xepr.acquire_scan()
plt.plot(run1.time,np.real(run1.data),label='run 1')
plt.plot(run2.time,np.real(run2.data),label='run 2')
plt.legend()
# plt.show()

t,nutation_data = exp.get_nutations(xepr,
                                    nu=[fc-0.3, fc+0.3],
                                    field=[fc,Bc],
                                    step=20e-3)
rp = resonatorProfile()
res_prof = rp.calc_res_prof(nutation_data, dt=t[1]-t[0])
rp.autofit([33,35])
rp.calculate_shape(res_prof,None)
rp.calc_IF_prof(2,rp.fc,type='fit')
fig = rp.res_prof_plot([33,35])
fig.savefig("resPro.png")
input("Press enter to move onto 5pDEER\n")

awg_exp.deer_pulse_5p(awg,100,[2.07,2.32],rp)

tau2 = input("5p Deer tau2?\n")


deadtime = 80
scans = 200
dim = np.floor((2*tau2 - deadtime)/16)


# Time to setup the Probe freq
xepr.set_freq(rp.fc)
print(f"Moving to the Probe Frequency:{rp.fc:.2f}\n")
exp.run_general(xepr,
    ["/PulseSpel/HUKA_DEER_AWG"],
    ["pi/2 tp=p1 +<x>/-<x>","+<x> (pi/2)"],
    {"PhaseCycle":False},
    {"p0":16,"p1":16,"h":20,"n":2})

xepr.run_exp()
input("Press Enter when the +<x> channel is correctly tuned?")

exp.run_general(xepr,
    ["/PulseSpel/HUKA_DEER_AWG"],
    ["pi/2 tp=p1 +<x>/-<x>","-<x> (pi/2)"],
    {"PhaseCycle":False},
    {"p0":16,"p1":16,"h":20,"n":2})

xepr.run_exp()
input("Press Enter when the -<x> channel is correctly tuned?")

print("Conducting a Field Sweep at probe frequency\n")
exp.run_general(xepr,
    ["/PulseSpel/HUKA_DEER_AWG"],
    ["Field Sweep +<x>/-<x> pg=200ns","Field sweep +<x>/-<x>"],
    {"PhaseCycle":False,"ReplaceMode":False},
    {"p0":16,"p1":16,"h":25,"n":1,"d0":600})

time.sleep(10)

fs_data = xepr.acquire_scan()
fs = FieldSweep()
fs.import_from_dataclass(fs_data)
fs.find_max()

xepr.set_field(fs.max_field)
Bc = xepr.get_field()
fc = xepr.get_counterfreq()
print(f"max field = {fs.max_field:.1f}")


# Time to setup the obs freq
obs_freq = rp.fc-0.07
xepr.set_freq(obs_freq)
print(f"Moving to the Obs Frequency:{obs_freq:.2f}\n")
print(f"Adjusting field tempoarily for tuning...\n")
xepr.set_field(Bc - 40)

exp.run_general(xepr,
    ["/PulseSpel/HUKA_DEER_AWG"],
    ["pi/2 tp=p1 +<x>/-<x>","+<x> (pi/2)"],
    {"PhaseCycle":False},
    {"p0":16,"p1":16,"h":20,"n":2})

xepr.run_exp()
input("Press Enter when the +<x> channel is correctly tuned?")

exp.run_general(xepr,
    ["/PulseSpel/HUKA_DEER_AWG"],
    ["pi/2 tp=p1 +<x>/-<x>","-<x> (pi/2)"],
    {"PhaseCycle":False},
    {"p0":16,"p1":16,"h":20,"n":2})

xepr.run_exp()
input("Press Enter when the -<x> channel is correctly tuned?")

exp.run_general(xepr,
    ["/PulseSpel/HUKA_DEER_AWG"],
    ["pi tp-p1 -<y>","-<y> echo up!"],
    {"PhaseCycle":False},
    {"p0":16,"p1":16,"h":20,"n":2})

xepr.run_exp()
print("Please tune for Positive Real (ECHO UP!)\n")
input("Press Enter when the -<x> channel is correctly tuned?")

print(f"Returning field to Bc\n")
xepr.set_field(Bc)


exp.run_general(xepr,
    ["/PulseSpel/HUKA_DEER_AWG"],
    ["5p DEER","DEER run AWG -+<x>"],
    {"PhaseCycle":True,"ReplaceMode":False},
    {"p0":16,"p1":16,"h":20,"n":scans,"d2":tau2,"d11":200,"d3":deadtime,"dim8":dim})

xepr.run_exp()
xepr.abort_exp()

input("Press check global phase\n")

awg.StartSigGen(3)
time.sleep(2)
xepr.run_exp()
