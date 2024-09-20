import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate as interp
from autodeer.colors import primary_colors



def calculate_optimal_tau(CPanalysis, MeasTime, SNR,target_step=0.015,target_shrt=None, ci=50, full_output=True):

    dataset = CPanalysis.dataset
    results = CPanalysis.fit_result
    averages = dataset.nAvgs * dataset.shots * dataset.nPcyc
    noise = results.noiselvl

    if target_shrt is None:
        target_shrt = dataset.attrs['reptime'] *1e-6 # us -> s
    nPointsInTime = lambda x: x * 3600 / target_shrt
    n_points = lambda x: x / (target_step*1e-3)
    dt=0.016

    x=np.linspace(1e-3, CPanalysis.axis.values.max(), 1000)
    y = np.linspace(0,0.2,1000)
    fit = results.evaluate(CPanalysis.fit_model, x)*results.scale
    # fitUncert = results.propagate(CPanalysis.fit_model, x)
    # fitUncertCi = fitUncert.ci(ci)*results.scale
    ub = CPanalysis.fit_model(x,*results.paramUncert.ci(ci)[:-1,1])*results.paramUncert.ci(ci)[-1,1]
    lb = CPanalysis.fit_model(x,*results.paramUncert.ci(ci)[:-1,0])*results.paramUncert.ci(ci)[-1,0]

    spl_fit_inverse = interp.InterpolatedUnivariateSpline(np.flip(fit/(noise*np.sqrt(averages)*np.sqrt(x*2/dt))),np.flip(x*2), k=3)
    spl_fit_inverse_lb = interp.InterpolatedUnivariateSpline(np.flip(lb/(noise*np.sqrt(averages)*np.sqrt(x*2/dt))),np.flip(x*2), k=3)
    spl_fit_inverse_ub = interp.InterpolatedUnivariateSpline(np.flip(ub/(noise*np.sqrt(averages)*np.sqrt(x*2/dt))),np.flip(x*2), k=3)

    optimal = spl_fit_inverse(SNR  / np.sqrt(nPointsInTime(MeasTime)))
    optimal_lb = spl_fit_inverse_lb(SNR / np.sqrt(nPointsInTime(MeasTime)))
    optimal_ub = spl_fit_inverse_ub(SNR/ np.sqrt(nPointsInTime(MeasTime)))

    if full_output:
        return optimal, optimal_lb, optimal_ub
    else:
        return optimal
    
def plot_optimal_tau(CPanlaysis, SNR,MeasTime=None, MaxMeasTime=24, labels=None, fig=None, axs=None):

    if fig is None and axs is None:
        fig, axs = plt.subplots(1,1)
    elif axs is None:
        axs = fig.add_subplot(111)
    
    MeasTimeAxis = np.logspace(0, np.log2(MaxMeasTime), 100, base=2)
    
    if not isinstance(SNR, (list, tuple)):
        SNR = [SNR]
    
    if labels is None:
        labels = [f'SNR = {snr}' for snr in SNR]

    for i,snr in enumerate(SNR):
        optimal, optimal_lb, optimal_ub = calculate_optimal_tau(CPanlaysis, MeasTimeAxis, snr, full_output=True)
        axs.plot(MeasTimeAxis, optimal, color=primary_colors[i], label=labels[i])
        axs.fill_between(MeasTimeAxis, optimal_lb, optimal_ub, color=primary_colors[i], alpha=0.3)

    axs.set_xlim(*axs.get_xlim())
    axs.set_ylim(*axs.get_ylim())

    if MeasTime is not None:
        if not isinstance(MeasTime, (list, tuple)):
            MeasTime = [MeasTime]
        if len(MeasTime) != len(SNR):
            raise ValueError('MeasTime and SNR must have the same length')
        for i, (mt, snr) in enumerate(zip(MeasTime, SNR)):
            optimal, optimal_lb, optimal_ub = calculate_optimal_tau(CPanlaysis, mt, snr, full_output=True)
            ylim = axs.get_ylim()
            axs.vlines(mt, *ylim,ls='--', color=primary_colors[i])
            xlim = axs.get_xlim()
            axs.hlines(optimal, *xlim,ls='--', color=primary_colors[i])
            axs.fill_between(xlim, (optimal_lb,optimal_lb), (optimal_ub,optimal_ub), color=primary_colors[i], alpha=0.3)
            


    axs.set_xlabel('Measurement time / h')
    axs.set_ylabel(r'$\tau_{evo}$ / $\mu$s')
    axs.legend(ncol=2, loc = 'lower right')

    return fig
