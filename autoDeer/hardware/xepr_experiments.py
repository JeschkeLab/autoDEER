import importlib
import time
import numpy as np
import re
import autoDeer.tools as tools
from scipy.optimize import minimize_scalar

MODULE_DIR = importlib.util.find_spec('autoDeer').submodule_search_locations[0]

# =============================================================================


def run_general(
        api, ps_file: tuple, exp: tuple, settings: dict, variables: dict,
        run: bool = True) -> None:
    """A function to run a general Pulse Spel experiment through autoDeer.

    Parameters
    ----------
    api : _type_
        The current Bruker Xepr class
    ps_file : tuple
        A tuple containing the file path to both the ".exp" and ".def" files.
    exp : tuple
        A tuple giving the name of the experiment and phase cycle.
    settings : dict
        A dictionary containing possible acquisition settings. Options include 
        ['ReplaceMode','PhaseCycle','Acquisition_mode']
    variables : dict
        A dictionary containg pulse spel variables to choose from can, these 
        can also be dimension of experiment.
    run : bool, optional
        Should the experiment run or just compile, by default True

    Raises
    ------
    ValueError
        If an input is of the wrong type.
    """
      
    if len(ps_file) == 1:
        # Assuming that both the EXP file and the DEF file have the same 
        # name bar-extention
        exp_file = MODULE_DIR + ps_file[0] + ".exp"
        def_file = MODULE_DIR + ps_file[0] + ".def"

    elif len(ps_file) == 2:
        
        # EXP and DEF file have seperate name
        exp_file = MODULE_DIR + ps_file[0] + ".exp"
        def_file = MODULE_DIR + ps_file[1] + ".def"

    else:
        raise ValueError(
            "ps_file must be of form ['EXP file'] or ['EXP file','DEF file']")

    # Identifying a dimension change in settings
    r = re.compile("dim([0-9]*)")
    match_list = list(filter(
        lambda list: list is not None, [r.match(i) for i in variables.keys()]))
    if len(match_list) >= 1:
        for i in range(0, len(match_list)):
            key = match_list[i][0]
            dim = int(r.findall(key)[0])
            new_length = int(variables[key])
            change_dimensions(exp_file, dim, new_length)
        
    api.set_PulseSpel_exp_filepath(exp_file)
    api.set_PulseSpel_def_filepath(def_file)
    api.compile_PulseSpel_prg()
    api.compile_PulseSpel_def()   

    if "ReplaceMode" in settings:
        api.set_ReplaceMode(settings["ReplaceMode"]) 
    else:
        api.set_ReplaceMode(False) 
    
    if "PhaseCycle" in settings:
        api.set_PhaseCycle(settings["PhaseCycle"]) 
    else:
        api.set_ReplaceMode(True) 

    if "Acquisition_mode" in settings:
        api.set_Acquisition_mode(settings["Acquisition_mode"])
    else:    
        api.set_Acquisition_mode(1)

    # setting PS Variables

    # Some Defaults first, these are overwritten if needed

    api.set_PulseSpel_var("p0", 16)
    api.set_PulseSpel_var("p1", 32)

    api.set_PulseSpel_var("d0", 400)
    api.set_PulseSpel_var("d1", 500)

    api.set_PulseSpel_var("d30", 16)
    api.set_PulseSpel_var("d31", 16)

    api.set_PulseSpel_var("h", 20)
    api.set_PulseSpel_var("n", 1000)
    api.set_PulseSpel_var("m", 1)

    # Change all variables

    for var in variables:
        api.set_PulseSpel_var(var.lower(), variables[var])

    api.set_PulseSpel_experiment(exp[0])
    api.set_PulseSpel_phase_cycling(exp[1])

    # Compile Defs and Program
    api.compile_PulseSpel_prg()
    api.compile_PulseSpel_def()  

    # Run Experiment
    if run is True:
        api.run_exp()
        time.sleep(1)
    pass

# =============================================================================


def change_dimensions(path, dim: int, new_length: int):    
    """A function to rewrite a pulseSpel experiment file with a new dimension

    Parameters
    ----------
    path : str
        The full file path.
    dim : int
        The experiment number that needs to be changed
    new_length : int
        The new length can be a list of two if 2D.

    Raises
    ------
    ValueError
        If there more than 2 dimesnions are supplied. Xepr can not handle 3+D 
        experiments.
    """
    # dim_str = "dim" + str(int(dim))
    # Open file
    with open(path, 'r') as file:
        data = file.readlines()

    # Build Regular expression to search for
    # This identifies the correct line and extracts the comment.
    re_search = fr"dim{int(dim)}\s*s\[[0-9]+,*[0-9]*\]" 
    if type(new_length) == int:
        new_string = f"dim{int(dim)} s[{int(new_length)}]"
    elif type(new_length) == list:
        if len(new_length) == 1:
            new_string = f"dim{int(dim)} s[{int(new_length[0])}]"
        elif len(new_length) == 2:
            new_string = \
                f"dim{int(dim)} s[{int(new_length[0])},{int(new_length[1])}]"
        else:
            raise ValueError("New length can't have more than 2 dimensions")
    else:
        raise ValueError("new_length must be either an int or a list")

    data = [re.sub(re_search, new_string, line) for line in data]

    with open(path, 'w') as file:
        file.writelines(data)

# =============================================================================


def get_nutations(api, nu, field, step, ELDOR: bool = True, nx: int = 128):

    min_freq = nu[0]
    max_freq = nu[1]

    freq_table = np.arange(min_freq, max_freq, step)

    n = len(freq_table)

    if len(field) == 2:
        input_freq = field[0]
        input_field = field[1]
        start_field = input_field * min_freq / input_freq
    elif len(field) == 1:
        start_field = field
    
    field_table = freq_table * start_field/min_freq
    
    # go to start field /  freq
    api.set_field(field_table[0], hold=True)
    api.set_freq(freq_table[0])
    if ELDOR:
        api.set_ELDOR_freq(freq_table[0])

    nut_data = np.zeros((n, nx+1), dtype=np.complex64)

    for i in range(0, n):
        api.set_field(field_table[i], hold=True)
        api.set_freq(freq_table[i])
        if ELDOR:
            api.set_ELDOR_freq(freq_table[i])

        api.run_exp()
        while api.is_exp_running():
            time.sleep(0.5)
        
        dataset = api.acquire_dataset()
        nut_data[i, 0] = api.get_counterfreq()
        nut_data[i, 1:] = dataset.data
        t = dataset.axes
        tools.progress_bar_frac(i, n)
    return t, nut_data

# =============================================================================


def CP_run(
        api, d0, num_pulses=3, ps_length=16, sweeps=4, dt=100, num_points=256):

    if num_pulses == 2:
        run_general(
            api,
            ["/PulseSpel/HUKA_DEER_AWG"],
            ["5p DEER relax", "DEER run AWG -+<x>"],
            {"PhaseCycle": True},
            {"p0": ps_length, "p1": ps_length, "h": 10, "n": sweeps, "d30": dt,
             "d0": d0, "dim10": num_points},
            False)
    else:
        raise ValueError("Only CP2 is currently implemented")

# =============================================================================


def DEER5p_run(
        api, ps_length, d0, tau2, sweeps=4, deadtime=80, dt=16,
        num_points=256):

    run_general(
        api,
        ["/PulseSpel/HUKA_DEER_AWG"],
        ["5p DEER", "DEER run AWG -+<x>"],
        {"PhaseCycle": True, "ReplaceMode": False},
        {"p0": ps_length, "p1": ps_length, "h": 20, "n": sweeps, "d2": tau2,
         "d11": 200, "d3": deadtime, "d30": dt, "d0": d0,
         "dim8": num_points},
        run=False)

# =============================================================================


class MPFUtune:
    """
    Tuning MPFU channels for optimal attenuation and phase
    """
    def __init__(self, api, echo="Hahn", ps_length=16, d0=680) -> None:
        """
        Parameters
        ----------
        api : _type_
            The spectrometr API object
        echo : str, optional
            The echo type. Options = ['Hahn","Refocused"], by default "Hahn"
        ps_length : int, optional
            The length of the pi/2 pulse, by default 16
        d0 : int, optional
            The approximate position of d0, this should be lower than ideal,
             by default 680
        """
        self.api = api
        self.hardware_wait = 5  # seconds 
        self.ps_length = ps_length
        self.d0 = d0
        if echo == "Hahn":
            self._setup_echo("Hahn Echo", tau1=400)
        elif echo == "Refocused":
            self._setup_echo("Refocused Echo", tau1=200, tau2=400)
        else:
            raise ValueError(
                "Only Hahn and Refocused echo's are currently supported")
        pass

    def _setup_echo(self, echo, tau1=400, tau2=400):
        PulseSpel_file = "/PulseSpel/phase_set"
        run_general(self.api,
                    [PulseSpel_file],
                    [echo, "BrXPhase"],
                    {"PhaseCycle": False},
                    {"p0": self.ps_length*2, "p1": self.ps_length, "h": 20,
                     "n": 1, "d0": self.d0, "d1": tau1, "d2": tau2, "pg": 128},
                    run=False
                    )

    def tune_phase(
            self, channel: str, target: str, tol=0.1, maxiter=30) -> float:
        """Tunes the phase of a given channel to a given target using the
        standard scipy optimisation scripts. 

        Parameters
        ----------
        channel : str
            The chosen MPFU channel. Options: ['+<x>', '-<x>', '+<y>', '-<y>']
        target : str
            The target echo position, this can either be maximising (+) or
            minimising (-) either the real (R) or imaginary (I) of the echo. 
            Options: ['R+', 'R-', 'I+', 'I-']
        tol : float, optional
            The tolerance in phase parameter, by default 0.1
        maxiter : int, optional
            The maximum number of iterations in the optimisation, by default 30

        Returns
        -------
        float
            The optimal value of the phase parameter

        """

        channel_opts = ['+<x>', '-<x>', '+<y>', '-<y>']
        phase_opts = ['R+', 'R-', 'I+', 'I-']

        if channel not in channel_opts:
            raise ValueError(f'Channel must be one of: {channel_opts}')
    
        if target not in phase_opts:
            raise ValueError(f'Phase target must be one of: {phase_opts}')
        if channel == '+<x>':
            phase_channel = 'BrXPhase'
        elif channel == '-<x>':
            phase_channel = 'BrMinXPhase'
        elif channel == '+<y>':
            phase_channel = 'BrYPhase'
        elif channel == '-<y>':
            phase_channel = 'BrMinYPhase'

        if target == 'R+':
            test_fun = lambda x: -1 * np.real(x)
        elif target == 'R-':
            test_fun = lambda x: 1 * np.real(x)
        elif target == 'I+':
            test_fun = lambda x: -1 * np.imag(x)
        elif target == 'I-':
            test_fun = lambda x: 1 * np.imag(x)

        lb = 0.0
        ub = 100.0

        def objective(x, *args):
            # x = x[0]
            self.api.hidden[phase_channel].value = x  # Set phase to value
            time.sleep(self.hardware_wait)
            self.api.run_exp()
            while self.api.is_exp_running():
                time.sleep(1)
            data = self.api.acquire_scan()
            v = data.data

            val = test_fun(np.sum(v))

            print(f'Phase Setting = {x:.1f} \t Echo Amplitude = {-1*val:.2f}')

            return val

        output = minimize_scalar(
            objective, method='bounded', bounds=[lb, ub],
            options={'xatol': tol, 'maxiter': maxiter})
        result = output.x
        print(f"Optimal Phase Setting for {phase_channel} is: {result:.1f}")
        self.api.hidden[phase_channel].value = result
        return result

    def tune_power(self, channel: str, tol=0.1, maxiter=30) -> float:
        """Tunes the attenuator of a given channel to a given target using the
        standard scipy optimisation scripts. 

        Parameters
        ----------
        channel : str
            The chosen MPFU channel. Options: ['+<x>', '-<x>', '+<y>', '-<y>']
        tol : float, optional
            The tolerance in attenuator parameter, by default 0.1
        maxiter : int, optional
            The maximum number of iterations in the optimisation, by default 30

        Returns
        -------
        float
            The optimal value of the attenuator parameter

        """
        channel_opts = ['+<x>', '-<x>', '+<y>', '-<y>']
        if channel not in channel_opts:
            raise ValueError(f'Channel must be one of: {channel_opts}')
        
        if channel == '+<x>':
            atten_channel = 'BrXAmp'
        elif channel == '-<x>':
            atten_channel = 'BrMinXAmp'
        elif channel == '+<y>':
            atten_channel = 'BrYAmp'
        elif channel == '-<y>':
            atten_channel = 'BrMinYAmp'

        lb = 0.0
        ub = 100.0

        def objective(x, *args):
            self.api.hidden[atten_channel].value = x  # Set phase to value
            time.sleep(self.hardware_wait)
            self.api.run_exp()
            while self.api.is_exp_running():
                time.sleep(1)
            data = self.api.acquire_scan()
            v = data.data

            val = -1 * np.sum(np.abs(v))

            print(f'Power Setting = {x:.1f} \t Echo Amplitude = {-1*val:.2f}')

            return val

        output = minimize_scalar(
            objective, method='bounded', bounds=[lb, ub],
            options={'xatol': tol, 'maxiter': maxiter})
        result = output.x
        print(f"Optimal Power Setting for {atten_channel} is: {result:.1f}")
        self.api.hidden[atten_channel].value = result
        return result

    def tune(self, channels: dict, tol: float = 0.1) -> None:
        """Tunes both the power and attenuation for a collection of channels.

        Parameters
        ----------
        channels : dict
            A dictionary of MPFU channels to be tunned and the associated phase
            target.\\
            Channel options = ['+<x>', '-<x>', '+<y>', '-<y>']\\
            Phase target options = ['R+', 'R-', 'I+', 'I-']\\
            E.g. {'+<x>': 'R+','-<x>': 'R-'}
        tol : float, optional
            The tolerance for all optimisations, by default 0.1
        """
        for channel in channels:
            if channel == '+<x>':
                phase_cycle = 'BrXPhase'
            elif channel == '-<x>':
                phase_cycle = 'BrMinXPhase'
            elif channel == '+<y>':
                phase_cycle = 'BrYPhase'
            elif channel == '-<y>':
                phase_cycle = 'BrMinYPhase'
            self.api.set_PulseSpel_phase_cycling(phase_cycle)
            print(f"Tuning channel: {channel}")
            self.tune_power(channel, tol=tol)
            self.tune_phase(channel, channels[channel], tol=tol)

    def calc_d0(self):
        initial_d0 = 500
        PulseSpel_file = "/PulseSpel/phase_set"
        run_general(
            self.api,
            [PulseSpel_file],
            ['Hahn Echo Trans', "BrXPhase"],
            {"PhaseCycle": False},
            {"p0": self.ps_length*2, "p1": self.ps_length, "h": 20, "n": 1,
             "d1": 400, "d0": initial_d0}
            )
        while self.api.is_exp_running():
            time.sleep(1)
        data = self.api.acquire_scan()
        max_pos = np.argmax(np.abs(data.data))
        max_time = data.axes[max_pos]
        d0 = initial_d0 + max_time
        return d0
