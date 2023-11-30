import numpy as np
import scipy.fft as fft
from scipy.interpolate import pchip_interpolate, splrep, BSpline
from scipy.signal import hilbert
import scipy.signal as signal
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pickle
import deerlab as dl


class ResonatorProfileAnalysis:

    def __init__(
            self, dataset,f_lims=(33,35)) -> None:
        """Analysis and calculation of resonator profiles.

        Parameters
        ----------
        dataset : xr.xarray
            The dataset containing the nutations. It must have both a 'LO' axis
            and a 'pulse0_tp' axis.
        f_lims : tuple, optional
            The frequency limits of the resonator profile, by default (33,35)
        """

        self.dataset = dataset.epr.correctphase

        self.freqs = self.dataset.LO
        self.n_files = self.freqs.shape[0]
        self.t = self.dataset.pulse0_tp
        self.f_lims = f_lims

        self._process_fit()
        pass

    def process_nutations(
            self, noisedensity: float = None, threshold:int = 2,
            nfft: int = 1000):
        """ Uses a power series to extract the resonator profile.

        Parameters
        ----------
        noisedensity : tuple, optional
            If not given the first trace is assumed to be so far off resonance
            that it is just noise. 
        nfft: int, optional
            The length of the fft to be used, zero padded if requred, default
            is 1000.
        threshold: int, optional
            The multiples above the noise a single must be to not be excluded,
            default is 2.

        Returns
        -------
        prof_data: np.ndarray
            The resonator profile, give in nutation frequency (GHz) 
        prof_frqs: np.ndarray
            The frequency axis in GHz

        """

        if noisedensity is None:
            f, Pxx_den = signal.welch(self.nutations[0, :], self.fs, nfft=nfft)

            noisedensity = np.mean(Pxx_den)
            
        threshold_val = threshold * noisedensity
        res_prof = np.zeros(self.n_files)
        for i in range(0, self.n_files):
            nuts = self.nutations[i, :]
            f, Pxx_den_sig = signal.welch(nuts, self.fs, nfft=nfft)
            nuoff = np.argwhere(Pxx_den_sig > threshold_val)
            if nuoff.shape[0] > 0:
                id = np.argmax(Pxx_den_sig[nuoff])
                res_prof[i] = np.abs(f[nuoff][id][0])
            else:
                res_prof[i] = 0
        
        nx_id = np.argwhere(res_prof > 0)
        self.prof_data = res_prof[nx_id]
        self.prof_frqs = np.real(self.freqs[nx_id])

        return self.prof_data, self.prof_frqs
    
    def _process_fit(self,R_limit=0.5):
        self.n_LO = self.freqs.shape[0]

        self.profile = np.zeros(self.n_LO)
        self.profile_ci = np.zeros(self.n_LO)

        fun = lambda x, f, tau,a,k: a*np.cos(2*np.pi*f*x)*np.exp(-x/tau)+ k
        bounds = ([5e-3,0,0,-1],[0.3,np.inf,2,1])

        R2 = lambda y, yhat: 1 - np.sum((y - yhat)**2) / np.sum((y - np.mean(y))**2)

        for i in range(self.n_LO):
            nutation = self.dataset[:,i]
            nutation = nutation/np.max(nutation)
            x = self.t
            try:
                results = curve_fit(fun, x, nutation, bounds=bounds,xtol=1e-4,ftol=1e-4)
                if R2(nutation,fun(x,*results[0])) > R_limit:
                    self.profile[i] = results[0][0]

                    self.profile_ci[i] = np.sqrt(np.diag(results[1]))[0]*1.96
                else:
                    self.profile[i] = np.nan
                    self.profile_ci[i] = np.nan
            except RuntimeError:
                self.profile[i] = np.nan
                self.profile_ci[i] = np.nan
            
        # Remove nan
        self.freqs = self.freqs[~np.isnan(self.profile)]
        self.profile = self.profile[~np.isnan(self.profile)]
        self.profile_ci = self.profile_ci[~np.isnan(self.profile_ci)]
        



    def fit(self,f_diff_threshold=2):
        """Fit the resonator profile with a sum of lorentzians.

        Parameters
        ----------
        f_diff_threshold : float, optional
            The difference between two peaks at which they will be merged into one, by default 0.03

        """

        frq_limits = {'lb': self.f_lims[0]-1, 'ub': self.f_lims[1]+1}
        def fit_gauss1():
            gauss_model = dl.dd_gauss
            gauss_model.mean.set(par0=34.05, **frq_limits)
            gauss_model.std.set(par0=0.2, lb=0.01, ub=10)

            result_gauss1 = dl.fit(gauss_model, self.profile, self.freqs,reg=False)
            return result_gauss1

        def fit_gauss2():
            gauss_model = dl.dd_gauss2
            gauss_model.mean1.set(par0=34.05, **frq_limits)
            gauss_model.std1.set(par0=0.2, lb=0.01, ub=10)
            gauss_model.mean2.set(par0=34.05, **frq_limits)
            gauss_model.std2.set(par0=0.3, lb=0.01, ub=10)

            result_gauss2 = dl.fit(gauss_model, self.profile, self.freqs,reg=False)
            return result_gauss2
            
        def fit_gauss3():
            gauss_model = dl.dd_gauss3
            gauss_model.mean1.set(par0=34.05, **frq_limits)
            gauss_model.std1.set(par0=0.2, lb=0.01, ub=10)
            gauss_model.mean2.set(par0=34.05, **frq_limits)
            gauss_model.std2.set(par0=0.3, lb=0.01, ub=10)
            gauss_model.mean3.set(par0=34.05, **frq_limits)
            gauss_model.std3.set(par0=0.3, lb=0.01, ub=10)

            result_gauss3 = dl.fit(gauss_model, self.profile, self.freqs,reg=False)
            return result_gauss3

        gs = [fit_gauss1(),fit_gauss2(),fit_gauss3()]
        n_modes = np.argmin([g.stats['chi2red'] for g in gs]) + 1

        print("Number of modes detected: ", n_modes)

        best_fit = gs[n_modes-1]
        # Find mean values which have an amplitude above 0.1 and are at least 0.1 apart
        # This is to avoid double counting
        if n_modes > 1:
            amps = np.array([getattr(best_fit,'amp'+str(i+1)) for i in range(n_modes)])
            means =  np.array([getattr(best_fit,'mean'+str(i+1)) for i in range(n_modes)])
            stds = np.array([getattr(best_fit,'std'+str(i+1)) for i in range(n_modes)])
        else:
            amps = np.array([best_fit.amp])
            means =  np.array([best_fit.mean])
            stds = np.array([best_fit.std])

        uni_means, uni_means_idx, uni_means_inverse = np.unique(np.around(means,decimals=f_diff_threshold),return_index=True,return_inverse=True)
        uni_amps = np.array([np.sum(amps[uni_means_inverse == i]) for i in np.unique(uni_means_inverse)])
        uni_stds = np.array([np.mean(stds[uni_means_inverse == i]) for i in np.unique(uni_means_inverse)])

        uni_means = uni_means[uni_amps > 0.01]
        uni_stds = uni_stds[uni_amps > 0.01]
        uni_amps = uni_amps[uni_amps > 0.01]

        # Remove modes which are outside the frequency range +- 0.5GHz
        uni_means = uni_means[(uni_means > self.f_lims[0]-0.5) & (uni_means < self.f_lims[1]+0.5)]
        uni_stds = uni_stds[(uni_means > self.f_lims[0]-0.5) & (uni_means < self.f_lims[1]+0.5)]
        uni_amps = uni_amps[(uni_means > self.f_lims[0]-0.5) & (uni_means < self.f_lims[1]+0.5)]

        n_modes = uni_means.shape[0]


        def lorenz_fcn(x, fc, q):
                    y = (0.5*fc)/(q*(x-fc)**2 + q*(0.5*(fc/q))**2)
                    return y

        def fit_lorenz1(gauss_fit,fit_kwargs):
            lorenz = dl.Model(lorenz_fcn, constants='x')
            lorenz.fc.par0 = uni_means[0]
            lorenz.fc.set(**frq_limits)
            lorenz.q.set(lb=0, ub=500)
            lorenz.q.par0 = uni_means[0]/uni_stds[0]
            lorenz.fc.set(par0=uni_means[0],lb=uni_means[0]-0.1,ub=uni_means[0]+0.1)
            result_lorenz1 = dl.fit(lorenz, self.profile, self.freqs,reg=False, **fit_kwargs)
            return result_lorenz1, lorenz_fcn, lorenz

        def fit_lorenz2(gauss_fit,fit_kwargs):
            def lorenz2_fcn(x, fc1, q1, fc2, q2, amp1, amp2):
                return lorenz_fcn(x, fc1, q1)*amp1 + lorenz_fcn(x, fc2, q2)*amp2

            lorenz2 = dl.Model(lorenz2_fcn,constants='x')
            lorenz2.fc1.par0 = uni_means[0]
            lorenz2.fc1.set(**frq_limits)
            lorenz2.q1.set(lb=0, ub=500)
            lorenz2.q1.par0 = uni_means[0]/uni_stds[0]
            lorenz2.fc2.par0 = uni_means[1]
            lorenz2.fc2.set(**frq_limits)
            lorenz2.q2.set(lb=0, ub=500)
            lorenz2.q2.par0 = uni_means[1]/uni_stds[1]
            lorenz2.amp1.set(lb=0, ub=np.inf, par0=uni_amps[0])
            lorenz2.amp2.set(lb=0, ub=np.inf, par0=uni_amps[1])
            result_lorenz2 = dl.fit(lorenz2, self.profile, self.freqs,reg=False, **fit_kwargs)
            return result_lorenz2, lorenz2_fcn, lorenz2

        def fit_lorenz3(gauss_fit,fit_kwargs):
            def lorenz3_fcn(x, fc1, q1, fc2, q2, fc3, q3, amp1, amp2, amp3):
                return lorenz_fcn(x, fc1, q1)*amp1 + lorenz_fcn(x, fc2, q2)*amp2 + lorenz_fcn(x, fc3, q3)*amp3
            lorenz3 = dl.Model(lorenz3_fcn,constants='x')
            lorenz3.fc1.par0 = uni_means[0]
            lorenz3.fc1.set(**frq_limits)
            lorenz3.q1.set(lb=0, ub=np.inf)
            lorenz3.q1.par0 = uni_means[0]/uni_stds[0]
            lorenz3.fc2.par0 = uni_means[1]
            lorenz3.fc2.set(**frq_limits)
            lorenz3.q2.set(lb=0, ub=np.inf)
            lorenz3.q2.par0 = uni_means[1]/uni_stds[1]
            lorenz3.fc3.par0 =  uni_means[2]
            lorenz3.fc3.set(**frq_limits)
            lorenz3.q3.set(lb=0, ub=np.inf)
            lorenz3.q3.par0 = uni_means[2]/uni_stds[2]
            lorenz3.amp1.set(lb=0, ub=np.inf, par0= uni_amps[0])
            lorenz3.amp2.set(lb=0, ub=np.inf, par0=uni_amps[1])
            lorenz3.amp3.set(lb=0, ub=np.inf, par0=uni_amps[2])

            result_lorenz3 = dl.fit(lorenz3, self.profile, self.freqs,reg=False, **fit_kwargs)
            return result_lorenz3, lorenz3_fcn, lorenz3

        self.model_x = np.linspace(self.f_lims[0], self.f_lims[1], 200)

        fit_kwargs = {'bootstrap': 100, 'bootcores':1}
        # fit_kwargs={}
        if n_modes == 1:
            lorenz_fit,fun,mod = fit_lorenz1(gs[0],fit_kwargs)
            self.model_func = lambda x: lorenz_fit.scale * lorenz_fcn(x, lorenz_fit.fc, lorenz_fit.q)
            self.fc = lorenz_fit.fc
            self.q = lorenz_fit.q
            
        elif n_modes == 2:
            lorenz_fit,fun,mod = fit_lorenz2(gs[1],fit_kwargs)
            self.model_func = lambda x: lorenz_fit.scale * fun(x, lorenz_fit.fc1, lorenz_fit.q1, lorenz_fit.fc2, lorenz_fit.q2, lorenz_fit.amp1, lorenz_fit.amp2)
            self.fc = np.array([lorenz_fit.fc1,lorenz_fit.fc2])
            self.q = np.array([lorenz_fit.q1,lorenz_fit.q2])

        elif n_modes == 3:
            lorenz_fit,fun,mod = fit_lorenz3(gs[2],fit_kwargs)
            self.model_func = lambda x: lorenz_fit.scale * fun(x, lorenz_fit.fc1, lorenz_fit.q1, lorenz_fit.fc2, lorenz_fit.q2, lorenz_fit.fc3, lorenz_fit.q3, lorenz_fit.amp1, lorenz_fit.amp2, lorenz_fit.amp3)
            self.fc = np.array([lorenz_fit.fc1,lorenz_fit.fc2,lorenz_fit.fc3])
            self.q = np.array([lorenz_fit.q1,lorenz_fit.q2,lorenz_fit.q3])

        self.model = self.model_func(self.model_x)
        self.lorenz_model = mod
        self.n_modes = n_modes
        self.pha = np.imag(hilbert(-np.log(self.model)))
        self.results = lorenz_fit

        return lorenz_fit

    def plot(self, fieldsweep=None, axs= None, fig=None):
        """plot. 

        Parameters
        ----------
        fieldsweep : FieldSweepAnalysis, optional
            Overlays the FieldSweep if provided, by default None
        axs : matplotlib.Axes, optional
            Axes to plot on, by default None
        fig : matplotlib.Figure, optional
            Figure to plot on, by default None

        Returns
        -------
        Matplotlib.Figure
            matplotlib figure object
        """
        if axs is None and fig is None:
            fig, axs = plt.subplots(2,1,constrained_layout=True,height_ratios=[0.8,0.2])

        if isinstance(axs, list) or isinstance(axs, np.ndarray):
            axs1 = axs[0]
            axs2 = axs[1]
        else:
            axs1 = axs
            axs2 = None
        prof_data = self.profile * 1e3 # GHz -MHz
        
        axs1.errorbar(self.freqs, prof_data, yerr=self.profile_ci*1e3, fmt='o', label='Data')

        if hasattr(self, 'model'):
            axs1.plot(self.model_x, self.model*1e3, label='Model')
            # CI= self.results.propagate(self.lorenz_model,self.model_x).ci(50)
            # axs1.fill_between(self.model_x, (CI[:,0])*1e3, (CI[:,1])*1e3, alpha=0.5)
            
            if axs2 is not None:
                axs2.plot(self.model_x, self.pha, label='Phase')
                axs2.set_xlabel("Frequency / GHz")
                axs2.set_ylabel("Phase / rad")
        


        axs1.set_xlabel("Frequency / GHz")
        axs1.set_ylabel("Nutation Frequency / MHz")

        def Hz2length(x):
            return 1 / ((x/1000)*2)

        def length2Hz(x):
            return 1 / ((x*1000)*2) 

        secax = axs1.secondary_yaxis('right', functions=(Hz2length, length2Hz))
        secax.set_ylabel(r'$\pi/2$ pulse length / ns')

        if fieldsweep:
            fsweep_data = np.abs(fieldsweep.data)
            fsweep_data /= fsweep_data.max()
            fsweep_data = fsweep_data * prof_data.max()
            axs1.plot(fieldsweep.fs_x + fieldsweep.LO, fsweep_data)

        return fig

# =============================================================================
# Spectral Position optimisation

def ceil(number, decimals=0):
    order = decimals*-1
    return np.ceil(number / 10**order)*10**order

def floor(number, decimals=0):
    order = decimals*-1
    return np.floor(number / 10**order)*10**order

def calc_overlap(x, func1, func2):
    """Calcuates the overlap between two functions.

    Parameters
    ----------
    x : np.ndarray
        The x axis of the functions
    func1 : function
        The first function
    func2 : function
        The second function

    Returns
    -------
    float
        The overlap between the two functions.
    """
    y1 = func1(x)
    y2 = func2(x)
    area_1 = np.trapz(y1, x)
    area_2 = np.trapz(y2, x)
    y1 /= area_1
    y2 /= area_2
    area_12 = np.trapz(y1*y2, x)
    return area_12

def BSpline_extra(tck_s):
    min_value = tck_s[0].min()
    max_value = tck_s[0].max()

    n_points = tck_s[1].shape[-1]
    n_points_10 = int(np.ceil(n_points/10))

    def pointwise(x):
        if x < min_value:
            return np.mean(tck_s[1][0:n_points_10])
        elif x > max_value:
            return np.mean(tck_s[1][-n_points_10:])
        else:
            return BSpline(*tck_s)(x)
        
    def ufunclike(xs):
        return np.array(list(map(pointwise, np.array(xs))))

    return ufunclike

def optimise_spectra_position(resonator_profile, fieldsweep, verbosity=0):

    band_width = resonator_profile.fc / resonator_profile.q

    if band_width > 0.4:
        n_points = 500
    else:
        n_points = 250

    upper_field_lim = ceil((resonator_profile.fc + band_width), 1)
    lower_field_lim = floor((resonator_profile.fc - band_width) , 1)
    shift_range = (upper_field_lim - lower_field_lim)/2

    x = np.linspace(lower_field_lim,upper_field_lim,n_points)
    # smooth_condition = (dl.der_snr(fieldsweep.data)**2)*fieldsweep.data.shape[-1]
    smooth_condition = 0.1
    x = np.flip(fieldsweep.fs_x + fieldsweep.LO)
    y = np.flip(fieldsweep.data)
    tck_s = splrep(x, y, s=smooth_condition)
    xc = signal.correlate(resonator_profile.fit_func(x),BSpline_extra(tck_s)(x), mode='same')
    x_shift = np.linspace(-1*shift_range,shift_range,n_points)
    new_LO = x_shift[xc.argmax()]

    return new_LO

