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


# =============================================================================
def DEERanalysis(dataset, compactness=True, model=None, ROI=False, verbosity=0, **kwargs):


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
        
        t = val_in_us(sequence.t)

    else:
        # Extract params from kwargs
        if "tau1" in kwargs:
            tau1 = kwargs["tau1"]
        if "tau2" in kwargs:
            tau2 = kwargs["tau2"]
        if "tau3" in kwargs:
            tau3 = kwargs["tau3"]
        exp_type = kwargs["exp_type"]
        t = dataset.axes[0]


    if verbosity > 1:
        print(f"")
        print(f"Experimental Parameters set to: tau1 {tau1:.2f} us \t tau2 {tau2:.2f} us")

    if "pathways" in kwargs:
        pathways = kwargs["pathways"]
    else:
        pathways = None
        
    if exp_type == "4pDEER":
        experimentInfo = dl.ex_4pdeer(tau1=tau1,tau2=tau2,pathways=pathways)
    elif exp_type == "5pDEER":
        experimentInfo = dl.ex_fwd5pdeer(tau1=tau1,tau2=tau2,tau3=tau3,pathways=pathways)
    elif exp_type == "3pDEER":
        experimentInfo = dl.ex_3pdeer(tau=tau1,pathways=pathways)

    r = np.linspace(1.5,10,100)


    if "regparamrange" in kwargs:
        regparamrange = kwargs["regparamrange"]
    else:
        regparamrange = None

    # identify experiment
    

    Vmodel = dl.dipolarmodel(t, r, experiment=experimentInfo)
    Vmodel.pathways = pathways

    if compactness:
        compactness_penalty = dl.dipolarpenalty(model, r, 'compactness', 'icc')
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

    if "mask" in kwargs:
        mask = kwargs["mask"]
        noiselvl = dl.noiselevel(Vexp[mask])
    else:
        noiselvl = dl.noiselevel(Vexp)
        mask=None


    # Core 
    if verbosity > 1:
        print('Starting Fitting')
    fit = dl.fit(Vmodel, Vexp, penalties=compactness_penalty, 
                 bootstrap=bootstrap, mask=mask, noiselvl=noiselvl,
                 regparamrange=regparamrange, bootcores=bootcores,verbose=verbosity)
    if verbosity > 1:
        print('Fit complete')
    mod_labels = re.findall(r"(lam\d*)'", str(fit.keys()))

    lam = 0
    for mod in mod_labels:
        lam += getattr(fit, mod)
    MNR = lam/fit.noiselvl
    fit.MNR = MNR
    fit.lam = lam


    fit.Vmodel = Vmodel
    fit.r = r
    fit.Vexp = Vexp
    fit.t = t
    fit.mask = mask
    
    if ROI:
        tau_max = lambda r: (r**3) *(2/4**3)
        fit.ROI = IdentifyROI(fit.P, r, criterion=0.90, method="gauss")
        rec_tau_max = tau_max(fit.ROI[1])
        return fit, rec_tau_max
    else:
        fit.ROI=None
        return fit
    


# =============================================================================

def DEERanalysis_plot(fit, background:bool, ROI=None, axs=None, fig=None):
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
def calc_overlap_frac(fieldsweep, pump_pulse, exc_pulse):
    
    f_axis,f_pump = pump_pulse._calc_fft(pad=1e4)
    f_pump = np.abs(f_pump)
    f_pump /= f_pump.max()
    
    _, f_exc = exc_pulse._calc_fft(pad=1e4)
    f_exc = np.abs(f_exc)
    f_exc /= f_exc.max()

    if np.iscomplexobj(fieldsweep.data):
        fs_data = correctphase(fieldsweep.data)
    else:
        fs_data = fieldsweep.data
    fs_data /= fs_data.max()

    fs_axis = fieldsweep.fs_x

    axis_min = np.fix(fieldsweep.fs_x.min()*1e2)/1e2
    axis_max = np.fix(fieldsweep.fs_x.max()*1e2)/1e2
    new_axis = np.linspace(axis_min, axis_max, 200)
    fs_new = interp(fs_axis, fs_data, new_axis)
    dead_spins = calc_overlap_region(
        f_axis, f_pump, f_axis, f_exc,
        new_axis=new_axis)

    total = np.trapz(fs_new, new_axis)
    num_pump = np.trapz(
        fs_new * interp(f_axis, f_pump, new_axis)
        - fs_new * dead_spins[1], new_axis)
    num_exc = np.trapz(
        fs_new * interp(f_axis, f_exc, new_axis) 
        - fs_new * dead_spins[1], new_axis)
    if np.isnan(num_exc):
        raise RuntimeError
    if np.isnan(num_pump):
        raise RuntimeError
    if total== 0:
        raise RuntimeError
    
    perc_pump = np.abs(num_pump)/total
    perc_exc = np.abs(num_exc)/total
    return calc_optimal_perc(perc_pump, perc_exc)



def calc_optimal_deer_frqs(
        fieldsweep, pump_pulse, exc_pulse, ref_pulse=None, pump_limits=None, exc_limits=None,
        method="density", ):
    """Calculates the necessary frequency shift to cause the maximal amount 
    of pump-exc overlap.

    Parameters
    ----------
    fieldsweep : FieldSweepAnalysis
        The fieldsweep of the spectrum. There is an assumption that the whole
        spectrum is spin 1/2
    pump_pulse : Pulse
        The pump pulse object.
    exc_pulse : Pulse
        The excitation pulse object.
    ref_pulse : Pulse
        The refocusing pulse object, this is used to adjust the efficacy of the 
        excitation pulse. This is of particular importance if the excitation pulse is
        wider than the refocusing pulse.
    pump_limits : tuple, optional
        The lower and upper limit of the pump frequency in GHz, by default 
        +- 0.3GHz
    exc_limits : tuple, optional
        The lower and upper limit of the pump frequency in GHz, by default 
        +- 0.3GHz
    method: str, optional
        Which method to calculate the excitation profiles. Either using
        a density matrix approach (slower, more accurate) or an FFT approach 
        (faster, less accurate). Options are "density" or "fft", by default
        "density"
    prof_res: float, optional
        The resolution of the excitation profile in MHz, by default 1.

    Returns
    -------
    pump frequency shift: float
        Frequency shift to make optimal of pump pulse in GHz.
    exc frequency shift: float
        Frequency shift to make optimal of excitation pulse in GHz.
    """

    if method == "fft":
        f_axis,f_pump = pump_pulse._calc_fft(pad=1e4)
        f_pump = np.abs(f_pump)
        f_pump /= f_pump.max()
        
        _, f_exc = exc_pulse._calc_fft(pad=1e4)
        f_exc = np.abs(f_exc)
        f_exc /= f_exc.max()
    elif method == "density":
        f_axis = np.arange(-1,1,0.005)
        _,_,f_pump = pump_pulse.exciteprofile(f_axis)
        f_pump = (f_pump - 1) * -0.5
        _,_,f_exc = exc_pulse.exciteprofile(f_axis)
        f_exc = (f_exc - 1) * -1
    
    if ref_pulse is not None:
        _,_,f_ref = ref_pulse.exciteprofile(f_axis)
        f_ref = (f_ref - 1) * -0.5
        f_exc *= f_ref
        f_exc *= f_ref

    if np.iscomplexobj(fieldsweep.data):
        fs_data = correctphase(fieldsweep.data)
    else:
        fs_data = fieldsweep.data
    fs_data /= fs_data.max()

    fs_axis = fieldsweep.fs_x

    def calc_optimal(pump_frq, exc_frq):
        axis_min = np.fix(fieldsweep.fs_x.min()*1e2)/1e2
        axis_max = np.fix(fieldsweep.fs_x.max()*1e2)/1e2
        new_axis = np.linspace(axis_min, axis_max, 200)
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
        if np.isnan(num_exc):
            raise RuntimeError
        if np.isnan(num_pump):
            raise RuntimeError
        if total== 0:
            raise RuntimeError
        
        perc_pump = np.abs(num_pump)/total
        perc_exc = np.abs(num_exc)/total
        return calc_optimal_perc(perc_pump, perc_exc)

    if pump_limits is not None:
        Pump_offset = np.linspace(pump_limits[0], pump_limits[1], 100)
    else:
        Pump_offset = np.linspace(-0.3, 0.3, 100)
    if exc_limits is not None:
        Exc_offset = np.linspace(exc_limits[0], exc_limits[1], 100)
    else:
        Exc_offset = np.linspace(-0.3, 0.3, 100)
    optimal_2d = np.array([[calc_optimal(p, e) for p in Pump_offset] for e in Exc_offset])

    # optimal_2d = np.abs(optimal_2d).reshape(100, 100)
    max_pos = np.unravel_index(np.nanargmax(optimal_2d), (100,100))
    print(
        f"Maximum achievable of {optimal_2d.max()*100:.2f}% optimal at"
        f"{Pump_offset[max_pos[1]]*1e3:.1f}MHz pump position and" 
        f" {Exc_offset[max_pos[0]]*1e3:.1f}MHz exc position.")

    return [Pump_offset[max_pos[1]], Exc_offset[max_pos[0]]]

def calc_optimal_perc(pump, exc):
    base_value = 0.5 * 0.5
    perc = (pump * exc) / base_value
    if np.isnan(perc):
        raise RuntimeError
    return perc 


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
        fieldsweep: FieldSweepAnalysis, pump_pulse,
        exc_pulse, respro: ResonatorProfileAnalysis = None,
        frq_shift: float = 0, axs=None, fig=None):
    """Generate a plot that contains the field sweep, pump and excitation 
    pulses as well as the resonator profile (optional). This should be used to 
    check that a DEER measurment is setup correctly and in an optimal
    configuration

    Parameters
    ----------
    fieldsweep : FieldSweepAnalysis
        The fieldsweep of the spectrum, with the frequency axis generated using
        the fieldsweep.calc_gyro(det_frq) command.
    pump_pulse : Pulse
        The pump pulse object.
    exc_pulse : int
        The excitation pulse object.
    respro : ResonatorProfileAnalysis, optional
        The resonator profile, by default None
    frq_shift : float, optional
        The frequency shift between resonator profile and fieldsweep maximum
        , by default 0
    
    """
    
    if fieldsweep.data.max() != 1.0:
        fieldsweep.data /= fieldsweep.data.max()

    # Build the pulses
    f_ax_pump,f_pump = pump_pulse._calc_fft(pad=1e4)
    f_pump = np.abs(f_pump)
    f_pump /= f_pump.max()
    
    f_ax_exc, f_exc = exc_pulse._calc_fft(pad=1e4)
    f_exc = np.abs(f_exc)
    f_exc /= f_exc.max()

    # Now build the plot
    axis_min = np.fix(fieldsweep.fs_x.min()*1e2)/1e2
    axis_max = np.fix(fieldsweep.fs_x.max()*1e2)/1e2
    new_axis = np.linspace(axis_min, axis_max, 200)
    fs_new = np.abs(interp(fieldsweep.fs_x, fieldsweep.data, new_axis))

    if axs is None and fig is None:
        fig, axs = plt.subplots()

    axs.plot(new_axis, fs_new, label="Field Sweep")
    dead_spins = calc_overlap_region(
        f_ax_pump, f_pump, f_ax_exc, f_exc, new_axis=new_axis)

    axs.fill(dead_spins[0], fs_new * dead_spins[1], color='0.6', label="Dead", 
            alpha=0.5, lw=0)
    axs.fill_between(
        new_axis, fs_new * interp(f_ax_exc, f_exc, new_axis),
        fs_new * dead_spins[1], label="Exc", color="C1", alpha=0.5, lw=0)
    axs.fill_between(
        new_axis, fs_new * interp(f_ax_pump, f_pump, new_axis),
        fs_new * dead_spins[1], label="Pump", color="C2", alpha=0.5, lw=0)
    if (respro is not None):
        axs.plot(respro.frq-respro.fc+frq_shift, respro.labs, label="Res Pro")
    axs.set_xlim(-0.25, 0.05)
    axs.legend()
    axs.set_ylabel(r"% flipped")
    axs.set_xlabel("Frequency / GHz")
    return fig

# =============================================================================
# Optimising pulse postions
# =============================================================================

def shift_pulse_freq(pulse,shift):
    if hasattr(pulse,'freq'):
        pulse.freq.value += shift
    elif hasattr(pulse,'BW'):
        if hasattr(pulse,'init_freq'):
            pulse.init_freq.value += shift
        elif hasattr(pulse,'final_freq'):
            pulse.final_freq.value += shift
    elif hasattr(pulse,'init_freq') and hasattr(pulse,'final_freq'):
        pulse.init_freq.value += shift
        pulse.final_freq.value += shift
    return pulse

def normalise_01(A):
    # Normalise the vector A to be between 0 and 1
    A = A - np.min(A)
    A = A / np.max(A)
    return A

def resample_and_shift_vector(A,f,shift):
    # Resample the vector A along axis f and shift it by shift and return on original axis f
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




def optimise_pulses(Fieldsweep, pump_pulse, exc_pulse, ref_pulse=None, filter=None, verbosity=0, method='grid',
                    nDEER=False, num_ref_pulses=2):
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

    Returns
    -------
    _type_
        _description_
    """


    gyro  = Fieldsweep.gyro
    fieldsweep_fun = lambda x: Fieldsweep.results.evaluate(Fieldsweep.model,(x+Fieldsweep.LO) /gyro*1e-1)

    f = np.linspace(-0.3,0.3,100)

    fieldsweep_profile = np.flip(fieldsweep_fun(f))
    
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
            fA,fB = np.unravel_index(Z.argmin(), Z.shape)
            if verbosity > 1:
                print(f"Filter: {filter:<8} Functional:{Z.min():.3f} \t Pump shift: {f[fA]*1e3:.1f}MHz \t Exc/Ref shift: {f[fB]*1e3:.1f}MHz")

        return Z.min(), f[fA], f[fB]
    
    if isinstance(filter, list):
        results = [optimise(f) for f in filter]
        best = np.argmin([r[0] for r in results])
        funct, fA, fB = results[best]
        if verbosity > 0:
            print(f"Best filter: {filter[best]}")
    elif filter is None:
        funct, fA, fB = optimise(None)
        if verbosity > 0:
            print(f"Filter: None Functional:{funct:.3f} \t Pump shift: {fA*1e3:.1f}MHz \t Exc/Ref shift: {fB*1e3:.1f}MHz")
    else:
        funct, fA, fB = optimise(filter)
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
        return new_pump_pulse, new_exc_pulse, new_ref_pulse, filter[best]
    else:
        return new_pump_pulse, new_exc_pulse, new_ref_pulse,


def plot_overlap(Fieldsweep, pump_pulse, exc_pulse, ref_pulse, filter=None, num_ref_pulses=2, axs=None, fig=None):
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
    num_ref_pulses : int, optional
        The total number of refocusing pulses, by default 2
    axs : matplotlib.axes, optional
        The axes to plot on, by default None
    fig : matplotlib.figure, optional
        The figure to plot on, by default None
    
    """


    gyro  = Fieldsweep.gyro
    fieldsweep_fun = lambda x: Fieldsweep.results.evaluate(Fieldsweep.model,(x+Fieldsweep.LO) /gyro*1e-1)

    f = np.linspace(-0.3,0.3,100)

    fieldsweep_profile = np.flip(fieldsweep_fun(f))

        
    pump_Mz = normalise_01(-1*pump_pulse.exciteprofile(freqs=f)[2].real)
    exc_Mz = normalise_01(-1*exc_pulse.exciteprofile(freqs=f)[2].real)

    if ref_pulse is not None:
        for i in range(num_ref_pulses):
            exc_Mz *= normalise_01(-1*ref_pulse.exciteprofile(freqs=f)[2].real)
    

    if filter == 'Matched':
        # Matched filter
        filter_profile = exc_Mz
    elif isinstance(filter, numbers.Number):
        filter_profile = build__lowpass_butter_filter(filter)

    if axs is None and fig is None:
        fig, axs = plt.subplots(1,1,figsize=(5,5))
    axs.plot(f,fieldsweep_profile, label = 'Fieldsweep', c='k')
    axs.fill_between(f,pump_Mz*fieldsweep_profile,linestyle='', label = 'Pump Profile', alpha=0.5,color='#D95B6F')
    axs.fill_between(f,exc_Mz*fieldsweep_profile, linestyle='', label = 'Observer Profile',alpha=0.5,color='#42A399')
    if filter is not None:
        axs.plot(f,filter_profile(f),'--', label = 'Filter')
    axs.legend()
    axs.set_xlabel('Frequency (MHz)')
    fig.tight_layout()

    return fig