import numpy as np
from scipy.optimize import curve_fit
import deerlab as dl
import matplotlib.pyplot as plt
from scipy.integrate import cumulative_trapezoid
import logging
import importlib
from autoDeer import FieldSweep
import scipy.fft as fft
from deerlab import correctphase
from scipy.interpolate import interp1d
import re


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
        compactness: bool = True, precision: str = "Detailed", **kwargs):
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

    if np.iscomplexobj(V):
        V = dl.correctphase(V)
    Vexp = V/np.max(V)         # Rescaling (aesthetic)
    t = t+zerotime

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
    
    if "pathways" in kwargs:
        pathways = kwargs["pathways"]
    if tau3 is not None:
        print("5pulse")
        experimentInfo = dl.ex_fwd5pdeer(tau1, tau2, tau3, pathways=pathways)
    else:
        experimentInfo = dl.ex_4pdeer(tau1, tau2, pathways=pathways)

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

    # Find total modulation depth
    mod_labels = re.findall("(lam\d*)'", str(fit.keys()))
    lam = 0
    for mod in mod_labels:
        lam += getattr(fit, mod)
    print(lam)
    MNR = lam/fit.noiselvl
    axs['Primary_time'].text(
        0.05, 0.05, f"MNR: {MNR:.2f}",
        transform=fig.transFigure, fontsize="16", color="black")
    if MNR < 10:
        axs['Primary_time'].text(
            0.15, 0.05, "LOW MNR: More averages requried",
            transform=fig.transFigure, fontsize="16", color="red")

    ROI_error = (rmax - ROI[1])
    
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

    return fig, ROI, rec_tau_max, fit

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


# =============================================================================

def calc_optimal_deer_frqs(
        fieldsweep: FieldSweep, pump_length: int, exc_length: int,
        det_frq: float = None, dt: float = 1):

    total_length = 1024
    pump_pulse = np.ones(int(pump_length/dt))
    pad_size = (total_length - pump_pulse.shape[0])/2
    pump_pulse = np.pad(
        pump_pulse, (int(np.floor(pad_size)), int(np.ceil(pad_size))))

    exc_pulse = np.ones(int(exc_length/dt))
    pad_size = (total_length - exc_pulse.shape[0])/2
    exc_pulse = np.pad(
        exc_pulse, (int(np.floor(pad_size)), int(np.ceil(pad_size))))

    f_axis = fft.fftshift(fft.fftfreq(total_length, dt))
    
    f_pump = np.abs(fft.fftshift(fft.fft(pump_pulse)))
    f_pump /= f_pump.max()
    
    f_exc = np.abs(fft.fftshift(fft.fft(exc_pulse)))
    f_exc /= f_exc.max()

    if np.iscomplexobj(fieldsweep.data):
        fs_data = correctphase(fieldsweep.data)
    fs_data /= fs_data.max()

    if det_frq is not None:
        fieldsweep.calc_gyro(det_frq)
    else:
        det_frq = fieldsweep.det_frq
    
    fs_axis = fieldsweep.fs_x

    def calc_overlap_region(x1, y1, x2, y2, new_axis=None):

        # Normalise data
        # y1 /= np.trapz(y1,x1)
        # y2 /= np.trapz(y2,x2)
        # create axis
        if new_axis is None:
            min_axis = np.max([np.min(x1), np.min(x2)])
            max_axis = np.min([np.max(x1), np.max(x2)])
            axis_mask = (x1 > min_axis) & (x1 < max_axis)
            x_new = x1[axis_mask]
            y1_new = y1[axis_mask]
        else:
            x_new = new_axis
            y1_interp = interp1d(x1, y1)
            y1_new = y1_interp(x_new)

        # Interpolate y2
        interp_y2 = interp1d(x2, y2)
        y2_new = interp_y2(x_new)

        return [x_new, np.minimum(y1_new, y2_new)]

    def interp(x1, y1, xnew):
        y1_interp = interp1d(x1, y1)
        y1_new = y1_interp(xnew)

        return y1_new

    def calc_optimal(pump_frq, exc_frq):
        new_axis = np.linspace(-0.27, 0.1, 200)
        fs_new = interp(fs_axis, fs_data, new_axis)
        dead_spins = calc_overlap_region(
            f_axis + pump_frq, f_pump, f_axis + exc_frq, f_exc,
            new_axis=new_axis)

        total = np.trapz(fs_new, new_axis)
        num_pump = np.trapz(
            fs_new * interp(f_axis + pump_frq, f_pump, new_axis)
            - fs_new * dead_spins[1], new_axis)
        num_exc = np.trapz(
            fs_new * interp(f_axis+exc_frq, f_exc, new_axis) 
            - fs_new * dead_spins[1], new_axis)

        perc_pump = num_pump/total
        perc_exc = num_exc/total
        return calc_optimal_perc(perc_pump, perc_exc)

    X = np.linspace(-0.15, 0.05, 100)
    Y = np.linspace(-0.15, 0.05, 100)
    optimal_2d = [calc_optimal(y, x) for y in X for x in Y]

    optimal_2d = np.abs(optimal_2d).reshape(100, 100)
    max_pos = np.unravel_index(optimal_2d.argmax(), (100, 100))
    print(
        f"Maxium achievable of {optimal_2d.max()*100:.2f}% optimal at"
        f"{X[max_pos[0]]*1e3:.1f}MHz pump position and" 
        f" {Y[max_pos[1]]*1e3:.1f}MHz exc position.")

    if det_frq is None:
        return [max_pos[0], max_pos[1]]
    else:
        return [X[max_pos[0]] + det_frq, Y[max_pos[1]] + det_frq] 


def calc_optimal_perc(pump, exc):
    base_value = 2/3 * np.sqrt(1/3)
    return (pump * np.sqrt(exc)) / base_value
