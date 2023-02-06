import numpy as np
from scipy.optimize import curve_fit
import deerlab as dl
import matplotlib.pyplot as plt
from scipy.integrate import cumulative_trapezoid
import logging
import importlib
from autodeer import FieldSweep, ResonatorProfile
import scipy.fft as fft
from deerlab import correctphase
from scipy.interpolate import interp1d
import re


log = logging.getLogger('core.DEER')


MODULE_DIR = importlib.util.find_spec('autodeer').submodule_search_locations

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
        compactness: bool = True, precision: str = "Detailed",
        background: bool = False, plot: bool = True, **kwargs):
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
    plot: bool, optional
        Should a figure be generated and returned, by default True.

    Returns
    -------
    _type_
        _description_
    """

    mask = None

    if np.iscomplexobj(V):
        V, Vim, ph = dl.correctphase(V, full_output=True)
        if "echos" in kwargs:
            echos = kwargs["echos"]
            mask = np.ones(V.shape, dtype=bool)
            for loc in echos:
                mask *= remove_echo(V, Vim, loc)
            
    if "mask" in kwargs:
        mask = kwargs["mask"]

    Vexp = V/np.max(V)         # Rescaling (aesthetic)
    t = t+zerotime

    if precision.lower() == "detailed":
        Ltype = 3
        regparamrange = [1e-8, 1e3]
    elif precision.lower() == "speed":
        Ltype = 6
        # regparamrange = [0.1, 1e3]
        regparamrange = [1e-2, 1e3]

    else:
        print("No Precision set. Ltype = 4nm")
        Ltype = 4
        regparamrange = [1e-8, 1e3]

    rmax = Ltype*(t.max()/2)**(1/3)
    rmax_e = np.ceil(rmax)  # Effective rmax
    
    r = np.linspace(1.5, rmax_e, num_points)  # nm  
    
    if "pathways" in kwargs:
        pathways = kwargs["pathways"]
    else:
        pathways = None
    if tau3 is not None:
        print("5pulse")
        pulses = 5
        experimentInfo = dl.ex_fwd5pdeer(tau1, tau2, tau3, pathways=pathways)
        if pathways is None:
            pathways = [1, 2, 3, 4, 5, 6, 7, 8]
    else:
        pulses = 4
        experimentInfo = dl.ex_4pdeer(tau1, tau2, pathways=pathways)
        if pathways is None:
            pathways = [1, 2, 3, 4]

    Vmodel = dl.dipolarmodel(t, r, experiment=experimentInfo)

    if (tau1 == tau2):
        # Link pathways
        if (pulses == 5) & (2 in pathways) & (8 in pathways):
            index2 = pathways.index(2) + 1
            index8 = pathways.index(8) + 1
            links = {"lam28": [f"lam{index2}", f"lam{index8}"]}
            Vmodel = dl.link(Vmodel, **links)
            pathways.remove(2)
            pathways.remove(8)
            pathways.append([2, 8])
        if (pulses == 4) & (1 in pathways) & (4 in pathways):
            index1 = pathways.index(1) + 1
            index4 = pathways.index(4) + 1
            Vmodel = dl.link(Vmodel, lam14=[f"lam{index1}", f"lam{index4}"])
            pathways.remove(1)
            pathways.remove(4)
            pathways.append([1, 4])

    Vmodel.pathways = pathways
    
    if compactness:
        compactness_penalty = dl.dipolarpenalty(None, r, 'compactness', 'icc')
        compactness_penalty.weight.set(lb=5e-3, ub=2e-1)
        # compactness_penalty.weight.freeze(0.04)

    else:
        compactness_penalty = None

    if "bootstrap" in kwargs:
        bootstrap = kwargs["bootstrap"]
    else:
        bootstrap = 0

    if "bootcores" in kwargs:
        bootcores = kwargs["bootcores"]
    else:
        bootcores = 1
    
    if mask is not None:
        noiselvl = dl.noiselevel(Vexp[mask])
    else:
        noiselvl = dl.noiselevel(Vexp)

    fit = dl.fit(Vmodel, Vexp, penalties=compactness_penalty, 
                 bootstrap=bootstrap, mask=mask, noiselvl=noiselvl,
                 regparamrange=regparamrange, bootcores=bootcores)

    mod_labels = re.findall(r"(lam\d*)'", str(fit.keys()))
    ROI = IdentifyROI(fit.P, r, criterion=0.90, method="int")

    fit.mask = mask

    # Find total modulation depth
    lam = 0
    for mod in mod_labels:
        lam += getattr(fit, mod)
    MNR = lam/fit.noiselvl
    fit.MNR = MNR
    
    rec_tau_max = (ROI[1]/3)**3 * 2

    # Expand fit object
    fit.Vmodel = Vmodel
    fit.r = r
    fit.Vexp = Vexp
    fit.t = t

    if plot:
        fig = DEER_analysis_plot(fit, background, ROI=ROI)
        return fit, ROI, rec_tau_max, fig
    else:
        return fit, ROI, rec_tau_max


# =============================================================================


def DEER_analysis_plot(fit, background, ROI=None):
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
    axs['Primary_time'].plot(t, Vfit, linewidth=3, color='C1', label='Fit')
    # axs['Primary_time'].fill_between(t, Vci[:, 0], Vci[:, 1], color='C1',
    #                                  alpha=0.3)
    if background:
        axs['Primary_time'].plot(
            t, background_func(t, fit), linewidth=3, color='C0',
            ls=':', label='Unmod. Cont.', alpha=0.5)

    axs['Primary_time'].set_xlabel(r"Time / $\mu s$")
    axs['Primary_time'].set_ylabel(r"A.U.")

    axs['Primary_time'].legend()

    # Distance Plots
    Pfit = fit.P
    Pci = fit.PUncert.ci(95)
    axs['Primary_dist'].plot(r, Pfit, '-', color='C1', lw=3, label='Fit')
    axs['Primary_dist'].fill_between(
        r, Pci[:, 0], Pci[:, 1], color='C1', alpha=0.3, label='95% CI')
    
    if ROI is not None:
        axs['Primary_dist'].fill_betweenx(
            [0, Pci[:, 1].max()], ROI[0], ROI[1], alpha=0.2, color='C2',
            label="ROI", hatch="/")
    
    axs['Primary_dist'].set_xlabel(r"Distance / $ nm$")
    axs['Primary_dist'].set_ylabel(r"$P(r^{-1})$")
    axs['Primary_dist'].legend()

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


# =============================================================================

def IdentifyROI(
        P: np.ndarray, r: np.ndarray, criterion: float = 0.99,
        method: str = "gauss"):
    """IdentifyROI _summary_

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
        Pmodel.mean.par0 = r[P.argmax()]
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
    else:
        fs_data = fieldsweep.data
    fs_data /= fs_data.max()

    if det_frq is not None:
        fieldsweep.calc_gyro(det_frq)
    else:
        det_frq = fieldsweep.det_frq
    
    fs_axis = fieldsweep.fs_x

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

    return [X[max_pos[0]], Y[max_pos[1]]]


def calc_optimal_perc(pump, exc):
    base_value = 2/3 * np.sqrt(1/3)
    return (pump * np.sqrt(exc)) / base_value


def interp(x1, y1, xnew):
    y1_interp = interp1d(x1, y1)
    y1_new = y1_interp(xnew)

    return y1_new


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


def plot_optimal_deer_frqs(
        fieldsweep: FieldSweep, pump_length: int, pump_frq: float,
        exc_length: int, exc_frq: float, respro: ResonatorProfile = None,
        frq_shift: float = 0, dt: float = 1):
    """Generate a plot that contains the field sweep, pump and excitation 
    pulses as well as the resonator profile (optional). This should be used to 
    check that a DEER measurment is setup correctly and in an optimal
    configuration

    Parameters
    ----------
    fieldsweep : FieldSweep
        The fieldsweep of the spectrum, with the frequency axis generated using
        the fieldsweep.calc_gyro(det_frq) command.
    pump_length : int
        The pump pulse length in ns.
    pump_frq : float
        The pump pulse frequnecy offset, from fieldsweep maximum in GHz
    exc_length : int
        The excitation pulse length in ns.
    exc_frq : float
        The excitation pulse frequnecy offset, from fieldsweep maximum in GHz
    respro : ResonatorProfile, optional
        The resonator profile, by default None
    frq_shift : float, optional
        The frequency shift between resonator profile and fieldsweep maximum
        , by default 0
    dt : float, optional
        The timestep in ns for generating pulses and their FFT, by default 1
    """
    
    if fieldsweep.data.max() != 1.0:
        fieldsweep.data /= fieldsweep.data.max()

    # Build the pulses
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

    # Now build the plot
    new_axis = np.linspace(-0.27, 0.1, 200)
    fs_new = interp(fieldsweep.fs_x, fieldsweep.data, new_axis)

    fig, ax = plt.subplots()
    ax.plot(new_axis, fs_new, label="Field Sweep")
    dead_spins = calc_overlap_region(
        f_axis+pump_frq, f_pump, f_axis+exc_frq, f_exc, new_axis=new_axis)

    ax.fill(dead_spins[0], fs_new * dead_spins[1], color='0.6', label="Dead", 
            alpha=0.5, lw=0)
    ax.fill_between(
        new_axis, fs_new * interp(f_axis+exc_frq, f_exc, new_axis),
        fs_new * dead_spins[1], label="Exc", color="C1", alpha=0.5, lw=0)
    ax.fill_between(
        new_axis, fs_new * interp(f_axis+pump_frq, f_pump, new_axis),
        fs_new * dead_spins[1], label="Pump", color="C2", alpha=0.5, lw=0)
    if (respro is not None):
        ax.plot(respro.frq-respro.fc+frq_shift, respro.labs, label="Res Pro")
    ax.set_xlim(-0.25, 0.05)
    ax.legend()
    ax.set_ylabel(r"% flipped")
    ax.set_xlabel("Frequency / GHz")