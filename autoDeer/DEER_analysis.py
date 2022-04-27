import time
import numpy as np
from scipy.optimize import curve_fit
import deerlab as dl
import matplotlib as plt
from autoDeer.DEER_4p import IdentifyROI
from mpl_toolkits.axes_grid1 import make_axes_locatable



def calc_identifiability(profile):
    x = profile['x']
    y = profile['y']
    fit_func = lambda x,a,b,c: a * (x-b)**2 + c
    fit_result = curve_fit(fit_func,x,y,bounds=(0,np.inf),xtol=0.001)
    return fit_result[0][0]


def std_4p_deer_analysis(t,V,tau1,tau2,zerotime=0,num_points=50,compactness=False):
    plt.style.use('seaborn')
    progress = 0
    progress_max = 6
    plt.ion()
    Vexp = dl.correctphase(V)
    Vexp = Vexp/np.max(Vexp)         # Rescaling (aesthetic)
    t = t-zerotime

    rmax = 6*(t.max()/2)**(1/3)
    rmax_e = np.ceil(rmax) # Effective rmax
    rmax_e = 5
    
    r = np.linspace(1,rmax_e,50) # nm  
    experimentInfo = dl.ex_4pdeer(tau1,tau2, pathways=[1,2,3])
    Vmodel = dl.dipolarmodel(t,r,experiment=experimentInfo)
    if compactness:
        compactness_penalty = dl.dipolarpenalty(None, r, 'compactness', 'icc')
    else:
        compactness_penalty = None
    fit = dl.fit(Vmodel,Vexp,penalties=compactness_penalty)
    ROI = IdentifyROI(fit.P,r,criterion=0.97)

    Vfit = fit.model
    Vci = fit.modelUncert.ci(95)

    fig, axs= plt.subplot_mosaic([
        ['Primary_time','Primary_time','Primary_dist','Primary_dist'],
        ['t_max_lam3','t_max_conc','r_max_lam3','r_max_conc']
    ],figsize=(15,8))
    fig.tight_layout(pad=4)
    fig.subplots_adjust(bottom=0.2,hspace=0.4)
    # Profile anaylsis
    Vmodel = dl.dipolarmodel(t,r,experiment=experimentInfo)
    num_t_profiles = int(np.ceil(t.max())- 1)
    end_points = np.around(-62.5*np.arange(0,num_t_profiles,1)).astype(int)
    end_points[0] = -1
    profiles_tmax = [None] * num_t_profiles
    for i in range(0,num_t_profiles):
        r = np.linspace(1.5,rmax_e,int(np.round(num_points/2))) # nm
        end_point = end_points[i]
        Vmodel = dl.dipolarmodel(t[0:end_point],r,experiment=experimentInfo)
        profiles_tmax[i] = dl.profile_analysis(Vmodel,Vexp[0:end_point],parameters=['lam3','conc'],
             grids={'lam3':np.linspace(0.2,0.7,10),'conc':np.linspace(50,300,10)})


    rmaxs = np.linspace(4,11,np.ceil((11-4)/2).astype(int)).astype(int)
    num_r_profiles = len(rmaxs)
    profiles_rmax = [None] * num_r_profiles

    for i in range(0,num_r_profiles):
        rmax = rmaxs[i]
        r = np.linspace(1.5,rmax,35)
        Vmodel = dl.dipolarmodel(t,r,experiment=experimentInfo)
        profiles_rmax[i] = dl.profile_analysis(Vmodel,Vexp,parameters=['lam3','conc'], 
            grids={'lam3':np.linspace(0.2,0.7,10),'conc':np.linspace(50,300,10)})

    # Top left Time domain
    axs['Primary_time'].plot(t,Vexp,'.',color='grey',label='Data')
    axs['Primary_time'].plot(t,Vfit,linewidth=3,color='orange',label='Fit')
    axs['Primary_time'].fill_between(t,Vci[:,0],Vci[:,1],color='orange',alpha=0.3)

    Bfcn = lambda mod,conc: (1-mod)*dl.bg_hom3d(time,conc,mod)
    Bfit = Bfcn(fit.lam3,fit.conc)
    axs['Primary_time'].plot(time,Bfit,'--',color='blue',label='Backgnd')


    axs['Primary_time'].set_xlabel(r"Time / $\mu s$")
    axs['Primary_time'].set_ylabel(r"A.U.")

    axs['Primary_time'].legend()
    # Top right distance domain
    Pfit = fit.P
    Pci = fit.PUncert.ci(95)
    r = np.linspace(1,rmax_e,50) # nm  
    axs['Primary_dist'].plot(r,Pfit,'-',color='orange',label='Fit')
    axs['Primary_dist'].fill_between(r,Pci[:,0],Pci[:,1],color='orange',alpha=0.3,label='95% CI')
    
    ROI = IdentifyROI(Pfit,r,criterion=0.97)
    axs['Primary_dist'].fill_betweenx([0,Pci[:,1].max()],ROI[0],ROI[1],alpha=0.2,color='green',label="ROI",hatch="/")
    
    axs['Primary_dist'].set_xlabel(r"Disatnce / $ nm$")
    axs['Primary_dist'].set_ylabel(r"$P(r^{-1})$")
    axs['Primary_dist'].legend()

    #Bottom Left Profile Analysis - Change Tmax

    colors = plt.cm.get_cmap('plasma')
    tst_lam3 = lambda x: calc_identifiability(x['lam3'].profile)
    tst_conc = lambda x: calc_identifiability(x['lam3'].profile)


    profiles_tmax_ident_lam3 = np.array(list(map(tst_lam3,profiles_tmax))) / 20
    profiles_rmax_ident_lam3 = np.array(list(map(tst_lam3,profiles_rmax))) / 20
    profiles_tmax_ident_conc = np.array(list(map(tst_conc,profiles_tmax))) / 20
    profiles_rmax_ident_conc = np.array(list(map(tst_conc,profiles_rmax))) / 20

    styles= ['dotted','dashed','dashdot','solid']
            

    for iD,profile in enumerate(profiles_tmax):
        profile = profile['lam3'].profile
        tmax = t[end_points[iD]]
        axs['t_max_lam3'].plot(profile['x'],profile['y'],label=fr"$t_{{max}}$ = {tmax:.1f} $\mu s$",
            color = colors(profiles_tmax_ident_lam3[iD]),linestyle = styles[iD])
 
    axs['t_max_lam3'].set_xlabel("lam3")
    axs['t_max_lam3'].set_ylabel("Profile Objective Function")
    axs['t_max_lam3'].set_title(r"$t_{max}$ profiles")
    axs['t_max_lam3'].legend(ncol=2,loc=(0,-0.4))

 

    #Bottom right Profile Analysis - Change rmax
    divider = make_axes_locatable(axs['r_max_conc'])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(plt.cm.ScalarMappable(cmap=colors), cax=cax)
    for iD,profile in enumerate(profiles_rmax):
        profile = profile['lam3'].profile
        rmax = rmaxs[iD]
        axs['r_max_lam3'].plot(profile['x'],profile['y'],label=fr"$r_{{max}}$ = {rmax:.1f} nm",
            color = colors(profiles_rmax_ident_lam3[iD]),linestyle = styles[iD])
 
    axs['r_max_lam3'].set_xlabel("lam3")
    axs['r_max_lam3'].set_ylabel("Profile Objective Function")
    axs['r_max_lam3'].legend(ncol=2,loc=(0,-0.4))
    axs['r_max_lam3'].set_title(r"$r_{max}$ profiles")
    
    #Bottom Left Profile Analysis - Change Tmax


            

    for iD,profile in enumerate(profiles_tmax):
        profile = profile['conc'].profile
        tmax = t[end_points[iD]]
        axs['t_max_conc'].plot(profile['x'],profile['y'],label=fr"$t_{{max}}$ = {tmax:.1f} $\mu s$",
        color = colors(profiles_tmax_ident_conc[iD]),linestyle = styles[iD])
 
    axs['t_max_conc'].set_xlabel("conc")
    axs['t_max_conc'].set_ylabel("Profile Objective Function")
    axs['t_max_conc'].set_title(r"$t_{max}$ profiles")
    axs['t_max_conc'].legend(ncol=2,loc=(0,-0.4))

 

    #Bottom right Profile Analysis - Change rmax
    for iD,profile in enumerate(profiles_rmax):
        profile = profile['conc'].profile
        rmax = rmaxs[iD]
        axs['r_max_conc'].plot(profile['x'],profile['y'],label=fr"$r_{{max}}$ = {rmax:.1f} nm",
            color = colors(profiles_rmax_ident_conc[iD]),linestyle = styles[iD])
 
    axs['r_max_conc'].set_xlabel("conc")
    axs['r_max_conc'].set_ylabel("Profile Objective Function")
    axs['r_max_conc'].legend(ncol=2,loc=(0,-0.4))
    axs['r_max_conc'].set_title(r"$r_{max}$ profiles")

    MNR = fit.lam3/fit.noiselvl
    axs['r_max_conc'].text(0.05,0.05,f"MNR: {fit.lam3/fit.noiselvl:.2f}",transform = fig.transFigure,size="x-large",color="black")
    if MNR < 20:
        axs['r_max_conc'].text(0.05,0.01,"LOW MNR: More averages requried",transform = fig.transFigure,size="x-large",color="red")

    ROI_error = (rmax - ROI[1])
    
    axs['r_max_conc'].text(0.55,0.05,f"ROI: {ROI[0]:.2f}nm to {ROI[1]:.2f}nm",transform = fig.transFigure,size="x-large",color="black")
    if ROI_error < 0.5:
        axs['r_max_conc'].text(0.55,0.01,r"ROI is too close to $r_{max}. Increase $\tau_{max}$",
        transform = fig.transFigure,size="x-large",color="red")



    fig.show()




    return [fit,profiles_tmax,profiles_rmax]

