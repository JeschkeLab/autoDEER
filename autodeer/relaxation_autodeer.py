import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate as interp
from autodeer.colors import primary_colors
from autodeer.Relaxation import CarrPurcellAnalysis, HahnEchoRelaxationAnalysis, RefocusedEcho2DAnalysis
from deerlab import der_snr


def calculate_optimal_tau(RelaxAnalysis, MeasTime, SNR, target_step=0.015, target_shrt=None, corr_factor=1, ci=50, full_output=True):
    """
    Calculate the optimal evolution time for a given SNR and measurement time.

    Parameters
    ----------
    RelaxAnalysis : CPanalysis
        The relaxation analysis object
    MeasTime : float or array-like
        The measurement time in hours. If array-like, the optimal tau will be calculated for each value.
    SNR : float
        The desired signal-to-noise ratio
    target_step : float, optional
        The target step size for the interpolation, by default 0.015 in microseconds
    target_shrt : float, optional
        The target shortest tau value for the interpolation in seconds. If None is given then the value is take from the relaxation experiment, by default None
    corr_factor : float, optional
        Correction factor for the MNR of the DEER fit result, by default 1
    ci : int, optional
        Confidence interval for fit in perc, by default 50
    full_output : bool, optional
        If True, the function will return the lower and upper bounds of the optimal tau, by default True
    
    Returns
    -------
    float or tuple
        The optimal tau value or a tuple with the lower and upper bounds of the optimal


    Notes
    -----
    * The uncertainty is not taken from deerlab as it must be a montomically decreasing function.

    """
    # Process the input data

    dataset = RelaxAnalysis.dataset
    averages = dataset.nAvgs * dataset.shots * dataset.nPcyc
    if target_shrt is None:
            target_shrt = dataset.attrs['reptime'] *1e-6 # us -> s

    if isinstance(RelaxAnalysis, RefocusedEcho2DAnalysis):
        # There will be no fit result, so we need to use the data
        if not hasattr(RelaxAnalysis, 'optimal_4p'):
            RelaxAnalysis._calc_optimal_4p_values()
        fit = RelaxAnalysis.optimal_4p_V
        axis = RelaxAnalysis.optimal_4p_tau2
        noise = der_snr(RelaxAnalysis.data/RelaxAnalysis.data.max())

        ub = fit + 2*noise
        lb = fit - 2*noise
        # enforce a lower bound on lb
        lb[lb<0] = 0

    else:

        if isinstance(RelaxAnalysis,CarrPurcellAnalysis):
            axis = RelaxAnalysis.axis.values * 2 # Tua_evo in us
        elif isinstance(RelaxAnalysis,HahnEchoRelaxationAnalysis):
            axis = RelaxAnalysis.axis.values
        else:
            raise ValueError('RelaxAnalysis must be an instance of CarrPurcellAnalysis, HahnEchoRelaxationAnalysis and RefocusedEcho2DAnalysis')
        data = RelaxAnalysis.data
        data /= data.max()
        R2 = 0
        if hasattr(RelaxAnalysis, 'fit_result'):
            fit_results = RelaxAnalysis.fit_result
            noise = fit_results.noiselvl

            # If the fit result is avaliable and off high quality, use it
            fit = fit_results.evaluate(RelaxAnalysis.fit_model, RelaxAnalysis.axis.values)*fit_results.scale
            R2 = fit_results.stats['R2']
            fitUncert = fit_results.propagate(RelaxAnalysis.fit_model, RelaxAnalysis.axis.values)
            # fitUncertCi = fitUncert.ci(50)*fit_results.scale
            # ub = fitUncertCi[:,1]
            # lb = fitUncertCi[:,0]
            ub = data + 2*noise
            lb = data - 2*noise
        if R2 < 0.9: # Either there is no fit or the fit is bad
            # Use the data
            data = RelaxAnalysis.data
            data /= data.max()
            noise = der_snr(data)
            fit = data
            ub = data + 2*noise
            lb = data - 2*noise
            # enforce a lower bound on lb
            lb[lb<0] = 0

    def find_all_numerical_roots(data,axis):
        sign_changes = np.where(np.abs(np.diff(np.sign(data)))>0)[0]
        if sign_changes.shape[0] == 0:
            if data[0] > 0:
                return axis[-1]
            else:
                return axis[0]
        else:
            return axis[sign_changes]

    nPointsInTime = lambda x: x * 3600 / target_shrt
    n_points = lambda x: x / (8*1e-3)


    functional = lambda SNR, T: corr_factor*fit - np.sqrt(SNR**2 * n_points(axis) * averages *noise**2 / (T*3600/target_shrt))
    functional_ub = lambda SNR, T: corr_factor*ub - np.sqrt(SNR**2 * n_points(axis) * averages *noise**2 / (T*3600/target_shrt))
    functional_lb = lambda SNR, T: corr_factor*lb - np.sqrt(SNR**2 * n_points(axis) * averages *noise**2 / (T*3600/target_shrt))

    if isinstance(MeasTime,np.ndarray):
        optimal = [find_all_numerical_roots(functional(SNR, T),axis).max() for T in MeasTime]
        optimal_lb = [find_all_numerical_roots(functional_lb(SNR, T),axis).max() for T in MeasTime]
        optimal_ub = [find_all_numerical_roots(functional_ub(SNR, T),axis).max() for T in MeasTime]
        optimal = np.array(optimal)
        optimal_lb = np.array(optimal_lb)
        optimal_ub = np.array(optimal_ub)
    else:
        optimal = find_all_numerical_roots(functional(SNR, MeasTime),axis).max()
        optimal_lb = find_all_numerical_roots(functional_lb(SNR, MeasTime),axis).max()
        optimal_ub = find_all_numerical_roots(functional_ub(SNR, MeasTime),axis).max()

    if full_output:
        return optimal, optimal_lb, optimal_ub
    else:
        return optimal
    
def plot_optimal_tau(CPanlaysis, SNR,MeasTime=None, MaxMeasTime=24, labels=None,corr_factor=1, fig=None, axs=None,cmap=None):
    """
    Generate a plot of the optimal evolution time for a given SNR against measurement time.

    Parameters
    ----------
    CPanlaysis : CPanalysis
        The analysis object for the Carr-Purcell Relaxation measurement
    SNR : float or array-like
        The desired signal-to-noise ratio. If array-like, a line will be drawn for each value.
    MeasTime : float or array-like, optional
        Draw crosshairs at specific measurement times, by default None. There must be the same number of elements in MeasTime and SNR
    MaxMeasTime : int, optional
        The maximum measurement time for the plot in hours, by default 24
    labels : list, optional
        Labels for the SNR values, by default None
    fig : matplotlib.figure.Figure, optional
        Figure to plot on, by default None
    axs : matplotlib.axes.Axes, optional
        Axes to plot on, by default None
    cmap : list, optional
        List of colors to use for the lines, by default None

    Returns
    -------
    matplotlib.figure.Figure
        The figure object
    
    """
    if fig is None and axs is None:
        fig, axs = plt.subplots(1,1)
    elif axs is None:
        axs = fig.add_subplot(111)
    
    if cmap is None:
        cmap = primary_colors

    MeasTimeAxis = np.logspace(0, np.log2(MaxMeasTime), 100, base=2)
    
    if not isinstance(SNR, (list, tuple)):
        SNR = [SNR]
    
    if labels is None:
        labels = [f'SNR = {snr}' for snr in SNR]

    for i,snr in enumerate(SNR):
        optimal, optimal_lb, optimal_ub = calculate_optimal_tau(CPanlaysis, MeasTimeAxis, snr, full_output=True, corr_factor=corr_factor)
        axs.plot(MeasTimeAxis, optimal, color=cmap[i], label=labels[i])
        axs.fill_between(MeasTimeAxis, optimal_lb, optimal_ub, color=cmap[i], alpha=0.3)

    axs.set_xlim(*axs.get_xlim())
    axs.set_ylim(*axs.get_ylim())

    if MeasTime is not None:
        if not isinstance(MeasTime, (list, tuple)):
            MeasTime = [MeasTime]
        if len(MeasTime) != len(SNR):
            raise ValueError('MeasTime and SNR must have the same length')
        for i, (mt, snr) in enumerate(zip(MeasTime, SNR)):
            optimal, optimal_lb, optimal_ub = calculate_optimal_tau(CPanlaysis, mt, snr, full_output=True,corr_factor=corr_factor)
            ylim = axs.get_ylim()
            axs.vlines(mt, *ylim,ls='--', color=cmap[i])
            xlim = axs.get_xlim()
            axs.hlines(optimal, *xlim,ls='--', color=cmap[i])
            axs.fill_between(xlim, (optimal_lb,optimal_lb), (optimal_ub,optimal_ub), color=cmap[i], alpha=0.3)
            


    axs.set_xlabel('Measurement time / h')
    axs.set_ylabel(r'$\tau_{evo}$ / $\mu$s')
    axs.legend(ncol=2, loc = 'lower right')

    return fig


def calc_correction_factor(CPanalysis,DEER_fit_result):
    """
    Calculate the correction factor for the MNR of the DEER fit result.

    Parameters
    ----------
    CPanalysis : CPanalysis
        CPanalysis object
    DEER_fit_result : DEERanalysis
        DEER_fit_result object
    
    Returns
    -------
    float
        Correction factor
    """
    shrt = DEER_fit_result.dataset.attrs['reptime'] *1e-6 # us -> s
    nPointsInTime = lambda x: x * 3600 / shrt
    t_axis = DEER_fit_result.t
    DEER_sequence = DEER_fit_result.dataset.epr.sequence

    actual_time = DEER_sequence._estimate_time() *(DEER_fit_result.dataset.nAvgs/DEER_fit_result.dataset.averages)

    # actual_MNR = DEER_fit_result.MNR
    actual_SNR = 1/DEER_fit_result.noiselvl

    if DEER_sequence.name == '5pDEER':
        dataset = CPanalysis.dataset
        CP_averages = dataset.nAvgs * dataset.shots * dataset.nPcyc
        noise = CPanalysis.fit_result.noiselvl
        V = lambda tau: CPanalysis.fit_result.evaluate(CPanalysis.fit_model, tau)[0]*CPanalysis.fit_result.scale

        tau = DEER_fit_result.dataset.attrs['tau1']/1e3

        est_SNR = V(tau)/(noise*np.sqrt(CP_averages))
        est_SNR *= np.sqrt(nPointsInTime(actual_time)/t_axis.shape[0])

        correction_factor = actual_SNR/est_SNR

    elif DEER_sequence.name == '4pDEER':
        dataset = CPanalysis.dataset
        CP_averages = dataset.nAvgs * dataset.shots * dataset.nPcyc
        noise = CPanalysis.fit_result.noiselvl
        V = lambda tau: CPanalysis.fit_result.evaluate(CPanalysis.fit_model, tau)[0]*CPanalysis.fit_result.scale

        tau = DEER_fit_result.dataset.attrs['tau2']/1e3

        est_SNR = V(tau)/(noise*np.sqrt(CP_averages))
        est_SNR *= np.sqrt(nPointsInTime(actual_time)/t_axis.shape[0])

        correction_factor = actual_SNR/est_SNR
    else:
        correction_factor = 1
    
    return correction_factor

