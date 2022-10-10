import numpy as np
from scipy.optimize import curve_fit
import deerlab as dl
import matplotlib.pyplot as plt
from scipy.integrate import cumulative_trapezoid
import logging
import importlib

log = logging.getLogger('core.DEER')


MODULE_DIR = importlib.util.find_spec('autoDeer').submodule_search_locations

# =============================================================================


def calc_identifiability(profile):
    x = profile['x']
    y = profile['y']
    fit_func = lambda x, a, b, c: a * (x-b)**2 + c
    fit_result = curve_fit(fit_func, x, y, bounds=(0, np.inf), xtol=0.001)
    return fit_result[0][0]

# =============================================================================


def std_deer_analysis(
        t: np.ndarray, V: np.ndarray, tau1: float, tau2: float,
        tau3: float = None, zerotime: float = 0, num_points: int = 50,
        compactness: bool = True, precision: str = "Detailed"):
    """std_deer_analysis This function conducts the standard deer analysis 
    using deerlab

    Parameters
    ----------
    t : np.ndarray
        The time domain of the data
    V : np.ndarray
        The signal, of equal length as t
    tau1 : float
        The first inter pulse delay.
    tau2 : float
        The second inter pulse delay
    tau3 : float, optional
        For 5p DEER the third delay. If this has a value it is automatically 
        detected as 5pulse, by default None
    zerotime : float, optional
        What value in the time axis represents the maximum, by default 0
    num_points : int, optional
        How many points in the distance distribution. This number should not 
        be too high, by default 50
    compactness : bool, optional
        Is the compactness criterion applied, this should always be applied, 
        by default True
    precision : str, optional
        How large the rmax should be. "Detailed", for normal work, "Speed" for 
        setup experiments, by default "Detailed"

    Returns
    -------
    _type_
        _description_
    """

    Vexp = dl.correctphase(V)
    Vexp = Vexp/np.max(Vexp)         # Rescaling (aesthetic)
    t = t-zerotime

    if precision == "Detailed":
        Ltype = 3
    elif precision == "Speed":
        Ltype = 6
    else:
        print("No Precision set. Ltype = 4nm")
        Ltype = 4

    rmax = Ltype*(t.max()/2)**(1/3)
    rmax_e = np.ceil(rmax)  # Effective rmax
    
    r = np.linspace(1, rmax_e, num_points)  # nm  
    
    experimentInfo = dl.ex_4pdeer(tau1, tau2, pathways=[1, 2, 3])
    
    if tau3 is not None:
        print("5pulse")
        experimentInfo = dl.ex_fwd5pdeer(tau1, tau2, tau3)
    else:
        experimentInfo = dl.ex_4pdeer(tau1, tau2, pathways=[1, 2, 3])

    Vmodel = dl.dipolarmodel(t, r, experiment=experimentInfo)
    
    if compactness:
        compactness_penalty = dl.dipolarpenalty(None, r, 'compactness', 'icc')
    else:
        compactness_penalty = None
    
    fit = dl.fit(Vmodel, Vexp, penalties=compactness_penalty)
    ROI = IdentifyROI(fit.P, r, criterion=0.99)

    Vfit = fit.model
    Vci = fit.modelUncert.ci(95)

    fig, axs = plt.subplot_mosaic([
        ['Primary_time', 'Primary_time', 'Primary_dist', 'Primary_dist']
        ], figsize=(12.5, 6.28))
    fig.tight_layout(pad=4)
    fig.subplots_adjust(bottom=0.2, hspace=0.4)

    # Top left Time domain
    axs['Primary_time'].plot(t, Vexp, '.', color='grey', label='Data')
    axs['Primary_time'].plot(t, Vfit, linewidth=3, color='orange', label='Fit')
    axs['Primary_time'].fill_between(
        t, Vci[:, 0], Vci[:, 1], color='orange', alpha=0.3)

    # Bfcn = lambda mod,conc: (1-mod)*dl.bg_hom3d(time,conc,mod)
    # Bfit = Bfcn(fit.lam3,fit.conc)
    # axs['Primary_time'].plot(time,Bfit,'--',color='blue',label='Backgnd')

    axs['Primary_time'].set_xlabel(r"Time / $\mu s$")
    axs['Primary_time'].set_ylabel(r"A.U.")

    axs['Primary_time'].legend()

    # Top right distance domain
    Pfit = fit.P
    Pci = fit.PUncert.ci(95)
    axs['Primary_dist'].plot(r, Pfit, '-', color='orange', label='Fit')
    axs['Primary_dist'].fill_between(
        r, Pci[:, 0], Pci[:, 1], color='orange', alpha=0.3, label='95% CI')
    
    axs['Primary_dist'].fill_betweenx(
        [0, Pci[:, 1].max()], ROI[0], ROI[1], alpha=0.2, color='green',
        label="ROI", hatch="/")
    
    axs['Primary_dist'].set_xlabel(r"Distance / $ nm$")
    axs['Primary_dist'].set_ylabel(r"$P(r^{-1})$")
    axs['Primary_dist'].legend()

    # **** Statistics ****
    MNR = fit.lam3/fit.noiselvl
    axs['Primary_time'].text(
        0.05, 0.05, f"MNR: {fit.lam3/fit.noiselvl:.2f}",
        transform=fig.transFigure, fontsize="16", color="black")
    if MNR < 10:
        axs['Primary_time'].text(
            0.05, 0.01, "LOW MNR: More averages requried",
            transform=fig.transFigure, fontsize="16", color="red")

    ROI_error = (rmax - ROI[1])
    
    rec_tau_max = (ROI[1]/3)**3 * 2

    axs['Primary_time'].text(
        0.55, 0.05, f"ROI: {ROI[0]:.2f}nm to {ROI[1]:.2f}nm",
        transform=fig.transFigure, fontsize="16", color="black")
    axs['Primary_time'].text(
        0.55, 0.01, rf"Recommended $\\tau_{{max}}$ = {rec_tau_max:.2f}us",
        transform=fig.transFigure, fontsize="16", color="black")

    if ROI_error < 0.5:
        axs['Primary_time'].text(
            0.55, 0.01,
            r"ROI is too close to $r_{max}. Increase $\\tau_{max}$",
            transform=fig.transFigure, size="x-large", color="red")

    fig.show()

    return fig, ROI, rec_tau_max

# =============================================================================


def IdentifyROI(
        dist: np.ndarray, r: np.ndarray, criterion: float = 0.99,
        plot: bool = False):
    """IdentifyROI _summary_

    Parameters
    ----------
    dist : np.ndarray
        The distance distribution.
    r : np.ndarray
        The distance axis
    criterion : float, optional
        The fraction of the distance distribution that must be in the ROI, by 
        default 0.99
    plot : bool, optional
        Plot the cumulative graphs, by default False
    """

    # Normlaize the distribution
    dist = dist / np.trapz(dist, r)
    # cumulative_dist = cumulative_trapezoid(dist,r,initial=0)
    # min_dist = r[np.argmin(np.abs(1 - cumulative_dist - criterion))]
    # max_dist = r[np.argmin(np.abs(cumulative_dist - criterion))]

    c_trapz_dist = np.zeros((dist.shape[0], dist.shape[0]))

    for i in range(0, dist.shape[0]):
        c_trapz_dist[i, i:] = cumulative_trapezoid(dist[i:], r[i:], initial=0)

    c_trapz_dist[(c_trapz_dist < criterion)] = 3
    ind = np.unravel_index(np.argmin(c_trapz_dist), c_trapz_dist.shape)
    min_dist = r[ind[0]]
    max_dist = r[ind[1]]

    # Enlarge ROI
    width = max_dist - min_dist
    max_dist = max_dist + width * 0.25
    min_dist = min_dist - width * 0.25

    if not plot:
        return [min_dist, max_dist]

    else:
        pass

# =============================================================================


def remove_echo(
        Vre: np.ndarray, Vim: np.ndarray, loc: int,
        Criteria: float = 4) -> np.ndarray:
    """This function removes crossing echoes. 
    Parameters
    ----------
    Vre : np.ndarray
        The real part of the phase corrected signal.
    Vim : np.ndarray
        The imaginary part of the phase corrected signal.
    loc : int
        The approximate location of the crossing echo, +- 20 data points
    Criteria : float, optional
        The delation criteria, in multiples of the std deviation, by default 4

    Returns
    -------
    np.ndarray
        The mask of points to be ignored.
    """

    search_mask = np.ones(Vre.shape[0], bool)
    search_mask[loc-30:loc+31] = False

    mask = np.abs(Vim) > Criteria * dl.noiselevel(Vre[search_mask])
    mask = mask & np.logical_not(search_mask)
    iskip = -1
    for i in range(len(mask)):
        if i < iskip:
            continue
        if mask[i]:
            mask[i-3:i] = True
            if i < len(mask) - 3:
                mask[i:i+3] = True
                iskip = i + 3

    mask = np.logical_not(mask)
    return mask
