import numpy as np
from scipy.optimize import curve_fit
import deerlab as dl
import matplotlib.pyplot as plt
from scipy.integrate import cumulative_trapezoid
import logging
import importlib
from autodeer import FieldSweepAnalysis, ResonatorProfileAnalysis
import scipy.fft as fft
from deerlab import correctphase
from scipy.interpolate import interp1d
import re
import numbers
import scipy.signal as sig
from scipy.optimize import minimize,brute
from autodeer.colors import primary_colors, ReIm_colors
from autodeer.utils import round_step


log = logging.getLogger('autoDEER.DEER')


MODULE_DIR = importlib.util.find_spec('autodeer').submodule_search_locations

# =============================================================================


def calc_identifiability(profile):
    x = profile['x']
    y = profile['y']
    fit_func = lambda x, a, b, c: a * (x-b)**2 + c
    fit_result = curve_fit(fit_func, x, y, bounds=(0, np.inf), xtol=0.001)
    return fit_result[0][0]

# =============================================================================

def find_longest_pulse(sequence):
    """
    Finds the longest pulse duration in a given sequence.
    
    Args:
    sequence (Sequence): The sequence to analyze.
    
    Returns:
    float: The duration of the longest pulse in microseconds.
    """
    longest_pulse = 0
    for pulse in sequence.pulses:
        tp = pulse.tp.value
        if tp > longest_pulse:
            longest_pulse = tp
    
    return longest_pulse/1e3

# =============================================================================
def DEERanalysis(dataset, compactness=True, model=None, ROI=False, exp_type='5pDEER', verbosity=0, **kwargs):


    Vexp:np.ndarray = dataset.data 

    if np.iscomplexobj(Vexp):
        Vexp = dl.correctphase(Vexp)

    Vexp /= Vexp.max()

    def val_in_us(Param):
        if len(Param.axis) == 0:
            if Param.unit == "us":
                return Param.value
            elif Param.unit == "ns":
                return Param.value / 1e3
        elif len(Param.axis) == 1:
            if Param.unit == "us":
                return sequence.tau1.value + Param.axis[0]['axis']
            elif Param.unit == "ns":
                return (Param.value + Param.axis[0]['axis']) / 1e3 

    if hasattr(dataset,"sequence"):
        sequence = dataset.sequence
        if sequence.name == "4pDEER":
            exp_type = "4pDEER"
            tau1 = val_in_us(sequence.tau1)
            tau2 = val_in_us(sequence.tau2)

        elif sequence.name == "5pDEER":
            exp_type = "5pDEER"
            tau1 = val_in_us(sequence.tau1)
            tau2 = val_in_us(sequence.tau2)
            tau3 = val_in_us(sequence.tau3)

        elif sequence.name == "3pDEER":
            exp_type = "3pDEER"
            tau1 = val_in_us(sequence.tau1)

        elif sequence.name == "nDEER-CP":
            exp_type = "4pDEER"
            tau1 = val_in_us(sequence.tau1)
            tau2 = val_in_us(sequence.tau2)

        tp = find_longest_pulse(sequence)
        t = val_in_us(sequence.t)
    
    elif 't' in dataset.coords: # If the dataset is a xarray and has sequence data
        t = dataset['t'].data
        if 'seq_name' in dataset.attrs:
            if dataset.attrs['seq_name'] == "4pDEER":
                exp_type = "4pDEER"
                tau1 = dataset.attrs['tau1'] / 1e3
                tau2 = dataset.attrs['tau2'] / 1e3
            elif dataset.attrs['seq_name'] == "5pDEER":
                exp_type = "5pDEER"
                tau1 = dataset.attrs['tau1'] / 1e3
                tau2 = dataset.attrs['tau2'] / 1e3
                tau3 = dataset.attrs['tau3'] / 1e3
            elif dataset.attrs['seq_name'] == "3pDEER":
                exp_type = "3pDEER"
                tau1 = dataset.attrs['tau1'] / 1e3
        else:
            if 'tau1' in dataset.attrs:
                tau1 = dataset.attrs['tau1'] / 1e3
            elif 'tau1' in kwargs:
                tau1 = kwargs.pop('tau1')
            else:
                raise ValueError("tau1 not found in dataset or kwargs")
            if 'tau2' in dataset.attrs:
                tau2 = dataset.attrs['tau2'] / 1e3
            elif 'tau2' in kwargs:
                tau2 = kwargs.pop('tau2')
            else:
                raise ValueError("tau2 not found in dataset or kwargs")
            if 'tau3' in dataset.attrs:
                tau3 = dataset.attrs['tau3'] / 1e3
                exp_type = "5pDEER"
            elif 'tau3' in kwargs:
                tau3 = kwargs.pop('tau3')
                exp_type = "5pDEER"
            else:
                exp_type = "4pDEER"

    else:
        # Extract params from kwargs
        if "tau1" in kwargs:
            tau1 = kwargs.pop("tau1")
        if "tau2" in kwargs:
            tau2 = kwargs.pop("tau2")
        if "tau3" in kwargs:
            tau3 = kwargs.pop("tau3")
        exp_type = kwargs.pop("exp_type")
        # t = dataset.axes[0]
        t = dataset['X']

    if exp_type == '4pDEER':
        pathways = [1,2,3]
    elif exp_type == '5pDEER':
        pathways = [1,5]
    elif exp_type == '3pDEER':
        pathways = [1,2]
    else:
        pathways = None

    if t.max() > 500:
        t /= 1e3

    if model is not None:
        compactness = False
    

    if verbosity > 1:
        print(f"Experiment type: {exp_type}, pathways: {pathways}")
        print(f"Experimental Parameters set to: tau1 {tau1:.2f} us \t tau2 {tau2:.2f} us")

    if "pathways" in kwargs:
        pathways = kwargs.pop("pathways")


    if hasattr(dataset,"sequence"):
        pulselength = tp
    elif hasattr(dataset,'attrs') and "pulse0_tp" in dataset.attrs:
        # seach over all pulses to find the longest pulse using regex
        pulselength = np.max([dataset.attrs[i] for i in dataset.attrs if re.match(r"pulse\d*_tp", i)])
        pulselength /= 1e3 # Convert to us
        pulselength /= 2 # Account for the too long range permitted by deerlab
    else:
        if "pulselength" in kwargs:
            pulselength = kwargs.pop("pulselength")/1e3
        else:
            pulselength = 16/1e3
        
    if exp_type == "4pDEER":
        experimentInfo = dl.ex_4pdeer(tau1=tau1,tau2=tau2,pathways=pathways,pulselength=pulselength)
    elif exp_type == "5pDEER":
        experimentInfo = dl.ex_fwd5pdeer(tau1=tau1,tau2=tau2,tau3=tau3,pathways=pathways,pulselength=pulselength)
    elif exp_type == "3pDEER":
        experimentInfo = dl.ex_3pdeer(tau=tau1,pathways=pathways,pulselength=pulselength)

  
    r_max = np.ceil(np.cbrt(t.max()*6**3/2))
    r = np.linspace(1.2,r_max,100)




    # identify experiment
    if 'parametrization' in kwargs:
        parametrization = kwargs.pop('parametrization')
    else:
        parametrization = 'reftimes'

    
    if 'Vmodel' in kwargs:
        Vmodel = kwargs['Vmodel']
    else:
        Vmodel = dl.dipolarmodel(t, r, experiment=experimentInfo, Pmodel=model,parametrization=parametrization)
        Vmodel.pathways = pathways

    if compactness:
        compactness_penalty = dl.dipolarpenalty(model, r, 'compactness', 'icc')
        compactness_penalty.weight.set(lb=5e-3, ub=2e-1)
        # compactness_penalty.weight.freeze(0.04)

    else:
        compactness_penalty = None

    if "mask" in kwargs:
        mask = kwargs["mask"]
        noiselvl = dl.noiselevel(Vexp[mask])
    else:
        noiselvl = dl.noiselevel(Vexp)
        mask=None


    # Cleanup extra args
    extra_args = ['tau1','tau2','tau3','exp_type','model','compactness','ROI','verbosity','pulselength','Vmodel','parametrization']
    for arg in extra_args:
        if arg in kwargs:
            kwargs.pop(arg)

    # Core 
    if verbosity > 1:
        print('Starting Fitting')
    fit = dl.fit(Vmodel, Vexp, penalties=compactness_penalty, 
                 noiselvl=noiselvl,
                 verbose=verbosity,
                 **kwargs)
    if verbosity > 1:
        print('Fit complete')
    mod_labels = re.findall(r"(lam\d*)'", str(fit.keys()))
    if mod_labels == []:
        mod_labels= ['mod']

    lam = 0
    for mod in mod_labels:
        lam += getattr(fit, mod)

    MNR = lam/fit.noiselvl
    fit.MNR = MNR
    fit.lam = lam


    fit.Vmodel = Vmodel
    fit.dataset = dataset
    fit.r = r
    fit.Vexp = Vexp
    fit.t = t
    fit.mask = mask

    if not hasattr(fit, "P"):
        fit.P = fit.evaluate(model,r)
        fit.PUncert = fit.propagate(model,r, lb=np.zeros_like(r))
    
    
    if ROI:
        # tau_max = lambda r: (r**3) *(3/4**3)
        tau_max = lambda r:  (r/3)**3 * 2
        fit.ROI = IdentifyROI(fit.P, r, criterion=0.90, method="gauss")
        if fit.ROI[0] < 1:
            fit.ROI[0] = 1
        rec_tau_max = tau_max(fit.ROI[1])
        fit.rec_tau_max = rec_tau_max
        dt_rec = lambda r: (r**3 *0.85)/(4*52.04)
        fit.rec_dt = np.max([dt_rec(fit.ROI[0]), 10e-3])

        return fit, rec_tau_max
    else:
        fit.ROI=None
        return fit
    


# =============================================================================
# Plotting
# =============================================================================
    
def background_func(t, fit):
    Vmodel = fit.Vmodel
    pathways = Vmodel.pathways
    conc = fit.conc
    prod = 1
    scale = 1
    if len(pathways) > 1:
        for i in pathways:
            if isinstance(i, list):
                # linked pathway
                lam = getattr(fit, f"lam{i[0]}{i[1]}")
                scale += -1 * lam
                for j in i:
                    reftime = getattr(fit, f"reftime{pathways.index(j)+1}")
                    prod *= dl.bg_hom3d(t-reftime, conc, lam)

            else:
                reftime = getattr(fit, f"reftime{i}")
                lam = getattr(fit, f"lam{i}")
                prod *= dl.bg_hom3d(t-reftime, conc, lam)
                scale += -1 * lam
    else:
        reftime = getattr(fit, f"reftime")
        lam = getattr(fit, f"mod")
        prod *= dl.bg_hom3d(t-reftime, conc, lam)
        scale += -1 * lam
    
    if hasattr(fit,'P_scale'):
        scale *= fit.P_scale
    elif hasattr(fit,'scale'):
        scale *= fit.scale

    return scale * prod

def calc_correction_factor(fit_result,aim_MNR=25,aim_time=2):
    """
    Calculate the correction factor for the number of averages required to achieve a given MNR in a given time.
    Parameters
    ----------
    fit_result : Deerlab.FitResult
        The fit result from the DEER analysis.
    aim_MNR : float, optional
        The desired MNR, by default 25
    aim_time : float, optional
        The desired time in hours, by default 2
    Returns
    -------
    float
        The correction factor for the number of averages.
    """

    dataset = fit_result.dataset
    runtime_s = dataset.nAvgs * dataset.nPcyc * dataset.shots * dataset.reptime * dataset.t.shape[0] * 1e-6
    aim_time *= 3600
    factor = fit_result.MNR /aim_MNR * np.sqrt(aim_time/runtime_s)
    return factor

def DEERanalysis_plot(fit, background:bool, ROI=None, axs=None, fig=None, text=True):
    """DEERanalysis_plot Generates a figure showing both the time domain and
    distance domain data along with extra important infomation such as the 
    Modulation to Noise Ratio (MNR), Region of Interest (ROI) and the 
    recommended dipolar evolution time for future experiments based upon the 
    ROI.

    Parameters
    ----------
    fit : Deerlab.FitResult
        _description_
    background : bool
        Should the background fit be plotted.
    ROI : tuple, optional
        The minimum and maximum of the Region of Interest (ROI),
        by default None

    Returns
    -------
    Figure
        A Matplotlib Figure object of the figure. 
    """
    mask = fit.mask
    t = fit.t
    r = fit.r
    Vexp = fit.Vexp
    Vfit = fit.model
    Vmodel = fit.Vmodel
    pathways = Vmodel.pathways

    # Calculate background
    def background_func(t, fit):
        conc = fit.conc
        prod = 1
        scale = 1
        for i in pathways:
            if type(i) is list:
                # linked pathway
                lam = getattr(fit, f"lam{i[0]}{i[1]}")
                scale += -1 * lam
                for j in i:
                    reftime = getattr(fit, f"reftime{pathways.index(j)+1}")
                    prod *= dl.bg_hom3d(t-reftime, conc, lam)

            else:
                reftime = getattr(fit, f"reftime{i}")
                lam = getattr(fit, f"lam{i}")
                prod *= dl.bg_hom3d(t-reftime, conc, lam)
                scale += -1 * lam

        return fit.P_scale * scale * prod

    if axs is None and fig is None:
        fig, axs = plt.subplot_mosaic([
            ['Primary_time', 'Primary_time', 'Primary_dist', 'Primary_dist']
            ], figsize=(12.5, 6.28))
        fig.tight_layout(pad=4)
        fig.subplots_adjust(bottom=0.2, hspace=0.4)

    # Time 

    if mask is not None:
        axs['Primary_time'].plot(
            t[mask], Vexp[mask], '.', color='0.7', label='Data', ms=6)
        axs['Primary_time'].plot(
            t[~mask], Vexp[~mask], '.', color='0.8', label='Data', ms=6)
    else:
        axs['Primary_time'].plot(t, Vexp, '.', color='0.7', label='Data', ms=6)
    axs['Primary_time'].plot(t, Vfit, linewidth=3, color=primary_colors[0], label='Fit')
    # axs['Primary_time'].fill_between(t, Vci[:, 0], Vci[:, 1], color='C1',
    #                                  alpha=0.3)
    if background:
        axs['Primary_time'].plot(
            t, background_func(t, fit), linewidth=3, color=primary_colors[1],
            ls=':', label='Unmod. Cont.', alpha=0.5)

    axs['Primary_time'].set_xlabel(r"Time / $\mu s$")
    axs['Primary_time'].set_ylabel(r"A.U.")

    axs['Primary_time'].legend()

    # Distance Plots
    Pfit = fit.P
    Pci = fit.PUncert.ci(95)
    axs['Primary_dist'].plot(r, Pfit, '-', color=primary_colors[0], lw=3, label='Fit')
    axs['Primary_dist'].fill_between(
        r, Pci[:, 0], Pci[:, 1], color=primary_colors[0], alpha=0.3, label='95% CI')
    
    if ROI is not None:
        axs['Primary_dist'].fill_betweenx(
            [0, Pci[:, 1].max()], ROI[0], ROI[1], alpha=0.2, color=primary_colors[1],
            label="ROI", hatch="/")
    
    axs['Primary_dist'].set_xlabel(r"Distance / $ nm$")
    axs['Primary_dist'].set_ylabel(r"$P(r) / nm^{-1}$")
    axs['Primary_dist'].legend()

    if not text:
        return fig
    # Analysis
    axs['Primary_time'].text(
        0.05, 0.05, f"MNR: {fit.MNR:.2f}",
        transform=fig.transFigure, fontsize="16", color="black")
    if fit.MNR < 10:
        axs['Primary_time'].text(
            0.15, 0.05, "LOW MNR: More averages requried",
            transform=fig.transFigure, fontsize="16", color="red")

    if ROI is not None:
        ROI_error = (r.max() - ROI[1])
    
        rec_tau_max = (ROI[1]/3)**3 * 2

        axs['Primary_time'].text(
            0.55, 0.05, f"ROI: {ROI[0]:.2f}nm to {ROI[1]:.2f}nm",
            transform=fig.transFigure, fontsize="16", color="black")
        axs['Primary_time'].text(
            0.05, 0.01, rf"Recommended $\tau_{{max}}$ = {rec_tau_max:.2f}us",
            transform=fig.transFigure, fontsize="16", color="black")

        if ROI_error < 0.5:
            axs['Primary_time'].text(
                0.55, 0.01,
                rf"ROI is too close to $r_{{max}}$. Increase $\tau_{{max}}$",
                transform=fig.transFigure, size="x-large", color="red")

    return fig

def DEERanalysis_plot_pub(results, ROI=None, fig=None, axs=None):
    """
    Generates a vertical plot of the DEER analysis results, ready for publication.
    
    Parameters
    ----------
    results : Deerlab.FitResult
        The results of the DEER analysis.
    ROI : tuple, optional
        The minimum and maximum of the Region of Interest (ROI),
        by default None
    fig : matplotlib.figure.Figure, optional
        The figure to plot the results on. If None, a new figure is created.
    axs : matplotlib.axes.Axes, optional
        The axes to plot the results on. If None, a new axes is created.
    """

    mask = results.mask
    t = results.t
    r = results.r
    Vexp = results.Vexp
    Vfit = results.model
    Vmodel = results.Vmodel
    pathways = Vmodel.pathways

    if fig is None:
        fig, axs = plt.subplots(2,1, figsize=(5,5),subplot_kw={},gridspec_kw={'hspace':0})
    else:
        axs = fig.subplots(2,1,subplot_kw={},gridspec_kw={'hspace':0})
    
    axs[0].plot(results.t,results.Vexp, '.',alpha=0.5,color='#D95B6F',mec='none')
    axs[0].plot(results.t,results.model, '-',alpha=1,color='#D95B6F', lw=2)
    
    axs[0].plot(results.t,background_func(results.t, results), '--',alpha=1,color='#42A399', lw=2)
    
    axs[0].set_xlabel('Time / $\mu s$')
    axs[0].xaxis.tick_top()

    axs[0].xaxis.set_label_position('top') 

    axs[0].set_ylabel('Signal A.U.')
    
    r = results.r
    axs[1].plot(r,results.P, '-',alpha=1.0,color='#D95B6F', label='P(r)')
    Pci = results.PUncert.ci(95)
    axs[1].fill_between(
        r, Pci[:, 0], Pci[:, 1], color='#D95B6F', alpha=0.3, label='P(r) 95% CI')
    if ROI is not None:
        axs[1].fill_betweenx(
            [0, Pci[:, 1].max()], ROI[0], ROI[1], alpha=0.2, color='#42A399',
            label="ROI", hatch="/")

    axs[1].set_xlabel('Distance / $nm$')
    axs[1].set_ylabel(r"$P(r) / nm^{-1}$")
    axs[1].legend()
    return fig

# =============================================================================

def IdentifyROI(
        P: np.ndarray, r: np.ndarray, criterion: float = 0.99,
        method: str = "gauss"):
    """IdentifyROI Identifies the region of interest. Two methods are sypported

    Methods
    +++++++
    
    1. Gaussian fitting ("gauss"):

    2. Intergration ("int"):

    Parameters
    ----------
    P : np.ndarray
        The distance distribution.
    r : np.ndarray
        The distance axis
    criterion : float, optional
        The fraction of the distance distribution that must be in the ROI, by 
        default 0.99
    method: str, optional
        The method used to calculate region of interest.
    """
    if method.lower() == "int":
        # Normlaize the distribution
        dist = P / np.trapz(P, r)
        # cumulative_dist = cumulative_trapezoid(dist,r,initial=0)
        # min_dist = r[np.argmin(np.abs(1 - cumulative_dist - criterion))]
        # max_dist = r[np.argmin(np.abs(cumulative_dist - criterion))]

        c_trapz_dist = np.zeros((dist.shape[0], dist.shape[0]))

        for i in range(0, dist.shape[0]):
            c_trapz_dist[i, i:] = cumulative_trapezoid(
                dist[i:], r[i:], initial=0)

        c_trapz_dist[(c_trapz_dist < criterion)] = 3
        ind = np.unravel_index(np.argmin(c_trapz_dist), c_trapz_dist.shape)
        min_dist = r[ind[0]]
        max_dist = r[ind[1]]

        # Enlarge ROI
        width = max_dist - min_dist
        max_dist = max_dist + width * 0.25
        min_dist = min_dist - width * 0.25

    elif method.lower() == "gauss":

        # Single gaussian approach
        Pmodel = dl.dd_gauss
        Pmodel.mean.lb=0
        Pmodel.mean.ub=r.max() *2
        Pmodel.mean.par0 = r[P.argmax()]
        Pmodel.std.lb=0
        Pmodel.std.ub=r.max()
        Pmodel.std.par0 = (r.max() - r.min())*0.5
        fit2 = dl.fit(Pmodel, P, r)
        min_dist = fit2.mean - 2.5 * fit2.std
        max_dist = fit2.mean + 2.5 * fit2.std

    return [min_dist, max_dist]


# =============================================================================

def remove_echo(
        Vre: np.ndarray, Vim: np.ndarray, loc: int,
        criteria: float = 4, extent: int = 3) -> np.ndarray:
    """This function removes crossing echoes. 
    Parameters
    ----------
    Vre : np.ndarray
        The real part of the phase corrected signal.
    Vim : np.ndarray
        The imaginary part of the phase corrected signal.
    loc : int
        The approximate location of the crossing echo, +- 30 data points
    criteria : float, optional
        The delation criteria, in multiples of the std deviation, by default 4
    extent : int, optional
        How many data points either side to remove, by default 3.

    Returns
    -------
    np.ndarray
        The mask of points to be ignored.
    """

    search_mask = np.ones(Vre.shape[0], bool)
    search_mask[loc-30:loc+31] = False

    mask = np.abs(Vim-np.mean(Vim)) > criteria * \
        dl.noiselevel(Vre[search_mask])
    mask = mask & np.logical_not(search_mask)
    iskip = -1
    for i in range(len(mask)):
        if i < iskip:
            continue
        if mask[i]:
            mask[i-extent:i] = True
            if i < len(mask) - extent:
                mask[i:i+extent] = True
                iskip = i + extent

    mask = np.logical_not(mask)
    return mask

# =============================================================================
# =============================================================================
# Optimising pulse postions
# =============================================================================

def shift_pulse_freq(pulse, shift):
    """
    Shifts the frequency of a pulse by a given amount.

    Args:
        pulse: The pulse whose frequency should be shifted.
        shift: The amount by which to shift the frequency.

    Returns:
        The pulse with the shifted frequency.
    """
    if hasattr(pulse, 'freq'):
        pulse.freq.value += shift
    elif hasattr(pulse, 'init_freq') and hasattr(pulse, 'final_freq'):
        pulse.init_freq.value += shift
        pulse.final_freq.value += shift
    return pulse

def normalise_01(A):
    """
    Normalizes the input vector A to be between 0 and 1.

    Parameters:
    A (numpy.ndarray): Input vector to be normalized.

    Returns:
    numpy.ndarray: Normalized vector between 0 and 1.
    """
    A = A - np.min(A)
    A = A / np.max(A)
    return A

def resample_and_shift_vector(A,f,shift):
    """
    Resample the vector A along axis f and shift it by shift and return on original axis f.

    Parameters:
    A (numpy.ndarray): The input vector to be resampled and shifted.
    f (numpy.ndarray): The axis along which to resample the vector.
    shift (float): The amount by which to shift the resampled vector.

    Returns:
    numpy.ndarray: The resampled and shifted vector.
    """
    A = np.interp(f,f+shift,A)
    return A


def build__lowpass_butter_filter(cutoff):
    """Build a lowpass butterworth filter with a cutoff frequency of cutoff

    Args:
        cutoff (float): cutoff frequency in GHz
    """

    sos = sig.butter(1,cutoff,btype='lowpass',analog=True)
    def filter_func(f):
        w = 2*np.pi*f
        w,h = sig.freqs(*sos,w)
        return abs(h)
    return filter_func

def functional(f_axis,fieldsweep,A,B,filter=None,A_shift=0,B_shift=0):
    """Functional for optimising the pulse positions

    Parameters
    ----------
    f_axis : np.ndarray
        The frequency axis of the field sweep in GHz
    fieldsweep : ad.FieldSweepAnalysis
        The FieldSweep analysis object
    A : np.ndarray
        The pump pulse profile
    B : np.ndarray
        The effective excitation pulse profile
    filter : np.ndarray, optional
        The filter profile if applicable, by default None
    A_shift : int, optional
        The shift in pump pulse in GHz, by default 0
    B_shift : int, optional
        The shift in effective exciatation pulse in GHz, by default 0

    Returns
    -------
    _type_
        _description_
    """
    D_profile = resample_and_shift_vector(A,f_axis,A_shift)*fieldsweep*resample_and_shift_vector(B,f_axis,B_shift)
    A_profile = fieldsweep*(resample_and_shift_vector(A,f_axis,A_shift))-D_profile # Pump
    B_profile = fieldsweep*(resample_and_shift_vector(B,f_axis,B_shift))-D_profile # Exc
    fsweep_norm = np.trapz(fieldsweep,f_axis)
    Aspins = np.trapz(A_profile,f_axis)/fsweep_norm

    if filter is None:
        Bspins = np.trapz(B_profile,f_axis)/fsweep_norm
        Noise = 1
    else:
        Bspins = np.trapz(B_profile*resample_and_shift_vector(filter,f_axis,B_shift),f_axis)/fsweep_norm
        Noise = np.trapz(filter,f_axis)

    return -1*(Aspins * Bspins)/np.sqrt(Noise)


def optimise_pulses(Fieldsweep, pump_pulse, exc_pulse, ref_pulse=None, filter=None, verbosity=0, method='brute',
                    nDEER=False, num_ref_pulses=2, full_output=False, resonator=None, **kwargs):
    """Optimise the pulse positions to maximise the pump-exc overlap.

    Parameters
    ----------
    Fieldsweep : ad.FieldSweepAnalysis
        The FieldSweep analysis object
    pump_pulse : ad.Pulse
        The pump pulse object
    exc_pulse : ad.Pulse
        The excitation pulse object
    ref_pulse : ad.Pulse, optional
        The refocusing pulse object\, by default None
    filter : str or number or list, optional
        The filter profile if applicable, by default None. If it is a number a filter is generated with this cutoff frequency.
        If the string 'Matched' is used a matched filter is used. If a list is used the optimisation is performed for each filter and the best is returned.
    verbosity : int, optional
        The verbosity, by default 0
    method : str, optional
        What search optimisation is used, by default 'grid'
    nDEER : bool, optional
        Is the sequence an nDEER sequrence, by default False. If True then the refocusing pulse is not optimised.
    num_ref_pulses : int, optional
        The total number of refocusing pulses, by default 2
    full_output : bool, optional
        Return the full output, by default False
    resonator : ad.ResonatorProfile, optional
        The resonator profile, by default None
    Returns
    -------
    ad.Pulse
        The optimised pump pulse
    ad.Pulse
        The optimised excitation pulse
    ad.Pulse
        The optimised refocusing pulse
    str or number
        The best filter, only if a list of filters is provided
    float
        The functional value after optimisation, only if full_output is True
    tuple
        The grid of the optimisation, only if full_output is True
    tuple
        The output of the optimisation, only if full_output is True
    """


    # gyro  = Fieldsweep.gyro
    # if hasattr(Fieldsweep,'results'):
    #     fieldsweep_fun = lambda x: Fieldsweep.results.evaluate(Fieldsweep.model,(x+Fieldsweep.LO) /gyro*1e-1)
    fieldsweep_fun = Fieldsweep.func_freq
    f = np.linspace(-0.3,0.3,100)
    fieldsweep_profile = fieldsweep_fun(f)

    # fieldsweep_profile = np.flip(fieldsweep_fun(f))
    
    pump_Mz = normalise_01(-1*pump_pulse.exciteprofile(freqs=f)[2].real)
    exc_Mz = normalise_01(-1*exc_pulse.exciteprofile(freqs=f)[2].real)

    if ref_pulse is not None:
        for i in range(num_ref_pulses):
            exc_Mz *= normalise_01(-1*ref_pulse.exciteprofile(freqs=f)[2].real)
    
    def optimise(filter):
        if filter == 'Matched':
            # Matched filter
            filter_profile = exc_Mz
        elif isinstance(filter, numbers.Number):
            filter_profile = build__lowpass_butter_filter(filter)(f)
        elif filter is None:
            filter_profile = None
            

        fun = lambda x: functional(f, fieldsweep_profile, pump_Mz, exc_Mz,filter_profile,A_shift=x[0],B_shift=x[1])
        
        if method == 'grid':
            X,Y = np.meshgrid(f,f)
            Z = [fun([x,y]) for x,y in zip(X.flatten(),Y.flatten())]
            Z = np.array(Z).reshape(X.shape)
            f_idA,f_idB = np.unravel_index(Z.argmin(), Z.shape)
            fA = f[f_idA]
            fB = f[f_idB]
            fval = Z.min()
            grid = (X,Y)
            Jout = Z
        elif method == 'brute':
            [fA,fB] ,fval,grid,Jout = brute(fun,(slice(f.min(),f.max(),0.01),slice(f.min(),f.max(),0.01)),full_output=True)
        
            if verbosity > 1:
                print(f"Filter: {filter:<8} Functional:{Z.min():.3f} \t Pump shift: {fA*1e3:.1f}MHz \t Exc/Ref shift: {fB*1e3:.1f}MHz")

        return fval, fA, fB,grid,Jout
    
    if isinstance(filter, list):
        results = [optimise(f) for f in filter]
        best = np.argmin([r[0] for r in results])
        funct, fA, fB = results[best]
        if verbosity > 0:
            print(f"Best filter: {filter[best]}")
    elif filter is None:
        funct, fA, fB,grid,Jout = optimise(None)
        if verbosity > 0:
            print(f"Filter: None Functional:{funct:.3f} \t Pump shift: {fA*1e3:.1f}MHz \t Exc/Ref shift: {fB*1e3:.1f}MHz")
    else:
        funct, fA, fB,grid,Jout  = optimise(filter)
        if verbosity > 0:
            print(f"Filter: {filter:<8} Functional:{funct:.3f} \t Pump shift: {fA*1e3:.1f}MHz \t Exc/Ref shift: {fB*1e3:.1f}MHz")

    new_pump_pulse = pump_pulse.copy()
    new_pump_pulse = shift_pulse_freq(new_pump_pulse,fA)
    new_exc_pulse = exc_pulse.copy()
    new_exc_pulse = shift_pulse_freq(new_exc_pulse,fB)

    new_ref_pulse = ref_pulse.copy()
    if not nDEER:
        new_ref_pulse = shift_pulse_freq(new_ref_pulse,fB)


    if isinstance(filter, list):
        if full_output:
            return new_pump_pulse, new_exc_pulse, new_ref_pulse, filter[best],funct,grid,Jout
        else:
            return new_pump_pulse, new_exc_pulse, new_ref_pulse, filter[best]
    else:
        if full_output:
            return new_pump_pulse, new_exc_pulse, new_ref_pulse, funct, grid, Jout
        else:
            return new_pump_pulse, new_exc_pulse, new_ref_pulse,


def plot_overlap(Fieldsweep, pump_pulse, exc_pulse, ref_pulse, filter=None, respro=None, num_ref_pulses=2, axs=None, fig=None):
    """Plots the pump and excitation profiles as well as the fieldsweep and filter profile.

    Parameters
    ----------

    Fieldsweep : ad.FieldSweepAnalysis
        The FieldSweep analysis object
    pump_pulse : ad.Pulse
        The pump pulse object
    exc_pulse : ad.Pulse
        The excitation pulse object
    ref_pulse : ad.Pulse, optional 
        The refocusing pulse object, by default None
    filter : str or number, optional
        The filter profile if applicable, by default None. If it is a number a filter is generated with this cutoff frequency.
        If the string 'Matched' is used a matched filter is used.
    respro : ad.ResonatorProfileAnalysis, optional
        The resonator profile for fitting, by default None. The resonator profile must include the fit.
    num_ref_pulses : int, optional
        The total number of refocusing pulses, by default 2
    axs : matplotlib.axes, optional
        The axes to plot on, by default None
    fig : matplotlib.figure, optional
        The figure to plot on, by default None
    
    """


    gyro  = Fieldsweep.gyro
    fieldsweep_fun = Fieldsweep.func_freq
    f = np.linspace(-0.4,0.4,100)

    fieldsweep_profile = fieldsweep_fun(f)

        
    pump_Mz = normalise_01(-1*pump_pulse.exciteprofile(freqs=f)[2].real)
    exc_Mz = normalise_01(-1*exc_pulse.exciteprofile(freqs=f)[2].real)

    if ref_pulse is not None:
        for i in range(num_ref_pulses):
            exc_Mz *= normalise_01(-1*ref_pulse.exciteprofile(freqs=f)[2].real)
    

    if filter == 'Matched':
        # Matched filter
        filter_profile  = lambda f_new: np.interp(f_new,f,exc_Mz)
    elif isinstance(filter, numbers.Number):
        filter_profile = build__lowpass_butter_filter(filter)

    if axs is None and fig is None:
        fig, axs = plt.subplots(1,1,figsize=(5,5), layout='constrained')
    elif axs is None:
        axs = fig.subplots(1,1,subplot_kw={},gridspec_kw={'hspace':0}
                           )
        
    # Normalise the fieldsweep profile
    fieldsweep_profile /= np.abs(fieldsweep_profile).max()
    
    # Plot the profiles
    axs.plot(f*1e3,fieldsweep_profile, label = 'Fieldsweep', c='k')
    axs.fill_between(f*1e3,pump_Mz*fieldsweep_profile, label = 'Pump Profile', alpha=0.5,color='#D95B6F')
    axs.fill_between(f*1e3,exc_Mz*fieldsweep_profile, label = 'Observer Profile',alpha=0.5,color='#42A399')
    if filter is not None:
        axs.plot(f*1e3,filter_profile(f),'--', label = 'Filter')

    if respro is not None:
        model_norm = respro.model / np.max(respro.model)
        axs.plot((respro.model_x - Fieldsweep.LO)*1e3, model_norm,'--', label='Resonator Profile')

    fmin = f[~np.isclose(fieldsweep_profile,0)].min()
    fmax = f[~np.isclose(fieldsweep_profile,0)].max()
    axs.set_xlim(fmin*1e3,fmax*1e3)

    axs.legend()
    axs.set_xlabel('Frequency (MHz)')

    return fig

def calc_deer_settings(experiment:str, CPdecay=None, Refocused2D=None, target_time=2,target_MNR=20, waveform_precision=2):
    """
    Calculates the optimal DEER settings based on the avaliable relaxation data

    Parameters
    ----------
    experiment : str
        Type of DEER experiment, either 'auto', '4pDEER' or '5pDEER'
    CPdecay : ad.CarrPurcellAnalysis
        Carr-Purcell relaxation data
    Refocused2D : ad.RefocusedEcho2DAnalysis, optional
        Refocused 2D data required for '4pDEER', by default None
    target_time : int, optional
        Target time for the DEER experiment in hours, by default 2
    target_MNR : float, optional
        Target modulation to noise ratio, by default 20
    waveform_precision : int, optional
        Precision of the waveform in ns, by default 2

    Returns
    -------
    dict
        DEER settings, with keys: 
            -'ExpType': '4pDEER' or '5pDEER'
            -'tau1': in us
            -'tau2': in us
            -'tau3': in us, only for 5pDEER
            -'AimTime': in hours

    Notes
    -----

    This function will calcate the optimal DEER settings based on the avaliable relaxation data, depending on the experiment type.
    For 4pDEER, the optimal tau1 and tau2 are calculated based on the refocused 2D data, and for 5pDEER, the optimal tau2 is calculated based on the CPdecay data or refocused 2D if CP decay data is not availiable.
    If the optimal tau2 for 5pDEER is less than 1.5us, the function will calculate the optimal tau1 and tau2 for 4pDEER instead. This is only possible if the refocused 2D data is availiable, otherwise a non optimal tau1 of 0.4us is used.

    """

    # Calculate the DEER settings from the avaliable relaxation data
    # Move out of the GUI Code asap

    deer_settings = {}

    if experiment == '4pDEER':
        if Refocused2D is None:
            raise ValueError("Refocused2D data required for 4pDEER")
        # Calcualting a 4pDEER Sequence
        deer_settings['ExpType'] = '4pDEER'

        # Calculate tau1 and tau2
        tau1,tau2 = Refocused2D.find_optimal(type='4pDEER',SNR_target=20, target_time=target_time, target_step=0.015)
        
        if tau2 < 2.5:
            # 2hr dipolar evo too short for meaningful results. Using double the time instead
            target_time *= 2
            tau1,tau2 = Refocused2D.find_optimal(type='4pDEER',SNR_target=20, target_time=target_time, target_step=0.015)
        
        if tau2 < 2.5:
            # Dipolar evo time still too short. Hardcoding a 2.5us dipolar evo time
            tau1,tau2 = Refocused2D.optimal_tau1(type='4pDEER',tau2=2.5)
        

        deer_settings['tau1'] = round_step(tau1,waveform_precision/1e3)
        deer_settings['tau2'] = round_step(tau2,waveform_precision/1e3)
        deer_settings['AimTime'] = target_time

    elif (experiment == '5pDEER') or (experiment == 'auto'):
        # Calculate a 5pDEER Sequence
        if CPdecay is not None:
            tau2hrs = CPdecay.find_optimal(SNR_target=target_MNR, target_time=target_time, target_step=0.015)
            tau4hrs = CPdecay.find_optimal(SNR_target=target_MNR, target_time=target_time*2, target_step=0.015)
        elif Refocused2D is not None:
            tau2hrs = Refocused2D.find_optimal(type='5pDEER',SNR_target=20, target_time=target_time, target_step=0.015)
            tau4hrs = Refocused2D.find_optimal(type='5pDEER',SNR_target=20, target_time=target_time*2, target_step=0.015)
        else:
            raise ValueError("CPdecay data required for 5pDEER")

        
        if (tau2hrs < 1.5) and (tau4hrs > 1.5):
            # 2hr dipolar evo too short for meaningful results. Using 4hr number instead
            deer_settings['tau1'] = round_step(tau4hrs,waveform_precision/1e3)
            deer_settings['tau2'] = round_step(tau4hrs,waveform_precision/1e3)
            deer_settings['tau3'] = round_step(0.3,waveform_precision/1e3)
            deer_settings['ExpType'] = '5pDEER'
            deer_settings['AimTime'] = target_time*2

        elif (tau2hrs < 1.5) and (tau4hrs < 1.5):
            # 2hr & 4hr dipolar evo too short for meaningful results. Hardcoding a 2.5us dipolar evo time, in 4pDEER and 4hrs
            # self.raise_warning(f"2hr dipolar evo too short. Hardcoding a 2.5us dipolar evo time")
            tau2 = 2.5
            if Refocused2D is not None:
                tau1 = Refocused2D.optimal_tau1(type='4pDEER',tau2=tau2)
            else:
                tau1 = 0.4

            deer_settings['tau1'] = round_step(tau1,waveform_precision/1e3)
            deer_settings['tau2'] = round_step(tau2,waveform_precision/1e3)
            
            deer_settings['ExpType'] = '4pDEER'
            deer_settings['AimTime'] = target_time * 2
            # self.raise_warning(f"2hr dipolar evo {tau2hrs:.2f} us, using 4pDEER")
        
        else:
            # Normal case, using 2hr dipolar evo time
            deer_settings['tau1'] = round_step(tau2hrs,waveform_precision/1e3)
            deer_settings['tau2'] = round_step(tau2hrs,waveform_precision/1e3)
            deer_settings['tau3'] = round_step(0.3,waveform_precision/1e3)
            deer_settings['ExpType'] = '5pDEER'
            deer_settings['AimTime'] = target_time

    return deer_settings
