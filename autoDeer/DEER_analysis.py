import time
import numpy as np
from scipy.optimize import curve_fit
import deerlab as dl
import matplotlib.pyplot as plt
from autoDeer.DEER_4p import IdentifyROI
from mpl_toolkits.axes_grid1 import make_axes_locatable



def calc_identifiability(profile):
    x = profile['x']
    y = profile['y']
    fit_func = lambda x,a,b,c: a * (x-b)**2 + c
    fit_result = curve_fit(fit_func,x,y,bounds=(0,np.inf),xtol=0.001)
    return fit_result[0][0]


def std_deer_analysis(t,V,tau1,tau2,tau3=None,zerotime=0,num_points=50,compactness=True):
    plt.style.use('seaborn')
    Vexp = dl.correctphase(V)
    Vexp = Vexp/np.max(Vexp)         # Rescaling (aesthetic)
    t = t-zerotime

    rmax = 6*(t.max()/2)**(1/3)
    rmax_e = np.ceil(rmax) # Effective rmax
    
    r = np.linspace(1,rmax_e,80) # nm  
    
    experimentInfo = dl.ex_4pdeer(tau1,tau2, pathways=[1,2,3])
    
    if tau3 != None:
        experimentInfo = dl.ex_fwd5pdeer(tau1,tau2,tau3)


    Vmodel = dl.dipolarmodel(t,r,experiment=experimentInfo)
    if compactness:
        compactness_penalty = dl.dipolarpenalty(None, r, 'compactness', 'icc')
    else:
        compactness_penalty = None
    fit = dl.fit(Vmodel,Vexp,penalties=compactness_penalty)
    ROI = IdentifyROI(fit.P,r,criterion=0.99)

    Vfit = fit.model
    Vci = fit.modelUncert.ci(95)

    fig, axs= plt.subplot_mosaic([
        ['Primary_time','Primary_time','Primary_dist','Primary_dist']
        ],figsize=(10,5))
    fig.tight_layout(pad=4)
    fig.subplots_adjust(bottom=0.2,hspace=0.4)

    # Top left Time domain
    axs['Primary_time'].plot(t,Vexp,'.',color='grey',label='Data')
    axs['Primary_time'].plot(t,Vfit,linewidth=3,color='orange',label='Fit')
    axs['Primary_time'].fill_between(t,Vci[:,0],Vci[:,1],color='orange',alpha=0.3)

    # Bfcn = lambda mod,conc: (1-mod)*dl.bg_hom3d(time,conc,mod)
    # Bfit = Bfcn(fit.lam3,fit.conc)
    # axs['Primary_time'].plot(time,Bfit,'--',color='blue',label='Backgnd')


    axs['Primary_time'].set_xlabel(r"Time / $\mu s$")
    axs['Primary_time'].set_ylabel(r"A.U.")

    axs['Primary_time'].legend()


    # Top right distance domain
    Pfit = fit.P
    Pci = fit.PUncert.ci(95)
    r = np.linspace(1,rmax_e,50) # nm  
    axs['Primary_dist'].plot(r,Pfit,'-',color='orange',label='Fit')
    axs['Primary_dist'].fill_between(r,Pci[:,0],Pci[:,1],color='orange',alpha=0.3,label='95% CI')
    
    axs['Primary_dist'].fill_betweenx([0,Pci[:,1].max()],ROI[0],ROI[1],alpha=0.2,color='green',label="ROI",hatch="/")
    
    axs['Primary_dist'].set_xlabel(r"Disatnce / $ nm$")
    axs['Primary_dist'].set_ylabel(r"$P(r^{-1})$")
    axs['Primary_dist'].legend()



    MNR = fit.lam3/fit.noiselvl
    axs['Primary_time'].text(0.05,0.05,f"MNR: {fit.lam3/fit.noiselvl:.2f}",transform = fig.transFigure,size="x-large",color="black")
    if MNR < 10:
        axs['Primary_time'].text(0.05,0.01,"LOW MNR: More averages requried",transform = fig.transFigure,size="x-large",color="red")

    ROI_error = (rmax - ROI[1])
    
    rec_tau_max = (ROI[1]/3)**3 * 2

    axs['Primary_time'].text(0.55,0.05,f"ROI: {ROI[0]:.2f}nm to {ROI[1]:.2f}nm",transform = fig.transFigure,size="x-large",color="black")
    axs['Primary_time'].text(0.55,0.01,rf"Recommended $\tau_{{max}}$ = {rec_tau_max:.2f}us",transform = fig.transFigure,size="x-large",color="black")

    if ROI_error < 0.5:
        axs['Primary_time'].text(0.55,0.01,r"ROI is too close to $r_{max}. Increase $\tau_{max}$",
        transform = fig.transFigure,size="x-large",color="red")

    fig.show()

    return fit

