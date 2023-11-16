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
            self, nuts: np.ndarray, freqs: np.ndarray, dt: float) -> None:
        """Analysis and calculation of resonator profiles.

        Parameters
        ----------
        nuts : np.ndarray
            An array containing the nutation data
        freqs: np.ndarray
            The respective frequencies of each nutation.
        dt : float
            The time step for the nutations.
        """
        if nuts.ndim != 2:
            raise ValueError(
                "nuts must be a 2D array where collumn 1 is the frequency")
        self.nutations = nuts
        self.freqs = freqs
        self.n_files = np.shape(nuts)[0]
        self.nx = np.shape(nuts)[1] - 1
        self.dt = dt
        self.fs = 1/dt
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

    def fit(self):
        def lorenz_fcn(x, centre, sigma):
            y = (0.5*sigma)/((x-centre)**2 + (0.5*sigma)**2)
            return y
        def lorenz_fcn(x, fc, q):
            y = (0.5*fc)/(q*(x-fc)**2 + q*(0.5*(fc/q))**2)
            return y


        self.prof_frqs = np.squeeze(self.prof_frqs)
        self.prof_data = np.squeeze(self.prof_data)

        lorenz = dl.Model(lorenz_fcn, constants='x')
        lorenz.fc.par0 = 34
        lorenz.fc.set(lb=33, ub=35)
        lorenz.q.set(lb=0, ub=100)
        lorenz.q.par0 = 80

        gauss_model = dl.dd_gauss
        gauss_model.mean.set(par0=34.05, lb=33, ub=35)
        gauss_model.std.set(par0=0.3, lb=0.01, ub=10)
        result_gauss = dl.fit(gauss_model, self.prof_data, self.prof_frqs)
        lorenz.fc.set(par0=result_gauss.mean,lb=result_gauss.mean-0.1,ub=result_gauss.mean+0.1)

        self.lorenz_model = lorenz
        results = dl.fit(lorenz, self.prof_data, self.prof_frqs)
        self.results = results

        self.fc = results.fc
        self.q = results.q

        self.profile_x = np.linspace(33, 35, 200)
        self.profile = lorenz_fcn(self.profile_x, self.fc, self.q)
        self.fit_func = lambda x: lorenz_fcn(x, self.fc, self.q)
        
        norm_prof = self.profile
        norm_prof /= norm_prof.max()
        self.pha = np.imag(hilbert(-np.log(self.profile)))

        pass

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
            fig, axs = plt.subplots(constrained_layout=True)
        prof_data = self.prof_data * 1e3 # GHz -MHz
        axs.plot(self.prof_frqs, prof_data)
        axs.set_xlabel("Frequency / GHz")
        axs.set_ylabel("Nutation Frequency / MHz")

        def Hz2length(x):
            return 1 / ((x/1000)*2)

        def length2Hz(x):
            return 1 / ((x*1000)*2) 

        secax = axs.secondary_yaxis('right', functions=(Hz2length, length2Hz))
        secax.set_ylabel(r'$\pi/2$ pulse length / ns')

        if fieldsweep:
            fsweep_data = np.abs(fieldsweep.data)
            fsweep_data /= fsweep_data.max()
            fsweep_data = fsweep_data * prof_data.max()
            axs.plot(fieldsweep.fs_x + fieldsweep.LO, fsweep_data)

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

