import autodeer.openepr as autoEPR
import numpy as np
import re
import time
import importlib
from autodeer import __version__
import os

MODULE_DIR = importlib.util.find_spec('autodeer').submodule_search_locations[0]

header = "; " + "-" * 79 + "\n"
header += "; " + \
     f"Auto-generated PulseSpel {{}} file by autoDEER {__version__}\n"
header += "; " + "-" * 79 + "\n"

# =============================================================================

class PSPhaseCycle:

    def __init__(self, sequence, MPFU=None, OnlyDet=False) -> None:
        BPhaseCycles = []
        
        if MPFU is not None:
            pulse_dicts = self._MPFU(sequence, MPFU)
        else:
            pulse_dicts = self._main(sequence)
        
        detect_dicts =self._detect(sequence)
        if not OnlyDet:
            BPhaseCycles += pulse_dicts
        BPhaseCycles += detect_dicts

        self.BPhaseCycles = BPhaseCycles

        self.__str__()

        pass

    def _MPFU(self, sequence, MPFU):
        BPhaseCycles = []
        for i,pulse in enumerate(sequence.pulses):
            if type(pulse) == autoEPR.Delay:
                continue
            elif type(pulse) == autoEPR.Detection:
                continue
            dict = {}
            dict["Pulse"] = i
            dict["Cycle"] = []
            for j in pulse.pcyc["Channels"]:
                dict["Cycle"].append(MPFU[j])
            BPhaseCycles.append(dict)
        
        return BPhaseCycles
    
    def _main(self, sequence):
        num_pulses = len(sequence.pcyc_vars)
        num_cycles = sequence.pcyc_dets.shape[0]

        BPhaseCycles = []
        # Normal pulses
        for i in range(0, num_pulses):
            dict = {}
            dict["Pulse"] = sequence.pcyc_vars[i]
            dict["Cycle"] = []
            for j in range(0, num_cycles):
                phase = sequence.pcyc_cycles[j, i]
                norm_phase = (phase / np.pi) % 2
                if norm_phase == 0:
                    B_phase = "+x"
                elif norm_phase == 1:
                    B_phase = "-x"
                elif norm_phase == 0.5:
                    B_phase = "+y"
                elif norm_phase == 1.5:
                    B_phase = "-y"
                dict["Cycle"].append(B_phase)
            BPhaseCycles.append(dict)

        return BPhaseCycles

    def _detect(self, sequence):
        num_cycles = sequence.pcyc_dets.shape[0]

        BPhaseCycles = []

        asg_dict = {"Cycle": [], "Pulse": 'a'}
        bsg_dict = {"Cycle": [], "Pulse": 'b'}
        for j in range(0, num_cycles):
            phase = sequence.pcyc_dets[j]
            if phase == 1:
                a_phase = "+a"
                b_phase = "+b"
            elif phase == -1:
                a_phase = "-a"
                b_phase = "-b"
            elif phase == 1j:
                a_phase = "+a"
                b_phase = "-b"
            elif phase == -1j:
                a_phase = "-a"
                b_phase = "+b"
            asg_dict["Cycle"].append(a_phase)
            bsg_dict["Cycle"].append(b_phase)

        BPhaseCycles.append(asg_dict)
        BPhaseCycles.append(bsg_dict)
        return BPhaseCycles



    def __str__(self) -> str:
        full_str = "begin lists \"auto\" \n"
        pcyc_hash = {}
        for i, cycle_dict in enumerate(self.BPhaseCycles):
            cycle = cycle_dict["Cycle"]
            if type(cycle_dict["Pulse"]) != str:
                pcyc_hash[cycle_dict["Pulse"]] = f"ph{i+1}"
                line_str = f"ph{i+1} "
            elif cycle_dict["Pulse"] == 'a':
                line_str = f"asg1 "
            elif cycle_dict["Pulse"] == 'b':
                line_str = f"bsg1 "
            for phase in cycle:
                line_str += f"{phase} "
            line_str += "\n"
            full_str += line_str

        full_str += "end lists \n"
        self.pcyc_hash = pcyc_hash
        return full_str


# =============================================================================


class PSparvar:

    def __init__(self, sequence: autoEPR.Sequence, id) -> None:

        progTable = sequence.progTable
        progTable_n = len(progTable["EventID"])

        parvar = {}
        self.PulseSpel = True
        parvar["variables"] = []
        self.events = []
        parvar["vec"] = []
        parvar["step"] = []

        for i in range(0, progTable_n):  # Loop over all progressions
            if progTable["axID"][i] == id:  # Ignore all progressions that aren't for this axis
                pulse_num = progTable["EventID"][i]
                var = progTable["Variable"][i]
                vec = progTable["axis"][i]
                parvar["variables"].append(var)
                parvar["vec"].append(vec)
                self.events.append(pulse_num)

                if var != "tp":
                    self.PulseSpel = False

                if len(np.unique(np.diff(vec))) != 1:
                    self.PulseSpel = False
                    parvar["step"].append(0)
                else:
                    parvar["step"].append(np.unique(np.diff(vec))[0])
        
        if "B" in parvar["variables"]:
            self.Bsweep = True
        else:
            self.Bsweep = False
        self.ax_step = parvar["step"][0]
        self.init = parvar["vec"][0][0]
        self.dim = parvar["vec"][0].shape[0]
        self.parvar = parvar
        pass

    def checkPulseSpel(self) -> bool:
        """Checks if parvar can be run in PulseSpel.\\
        Criteria:
            1) Only pulse/delay lengths changing
            2) Constant integer step changes

        Returns
        -------
        bool
            _description_
        """


possible_delays = [
    "d1", "d2", "d3", "d4", "d5", "d6", "d7", "d8", "d10", "d11", "d12",
    "d13", "d14", "d15", "d16", "d17", "d18", "d19", "d20", "d21", "d22",
    "d23", "d24", "d25", "d26", "d27", "d28", "d29", "d30"]

possible_pulses = [
        "p0", "p1", "p2", "p3", "p4", "p5", "p6", "p7", "p8", "p9", "p10"]

# =============================================================================


class PulseSpel:

    def __init__(self, sequence, MPFU=None, AWG=False) -> None:
        self.possible_delays = possible_delays.copy()
        self.possible_pulses = possible_pulses.copy()
        self.sequence = sequence
        self.var_hash = {}
        self.def_file_str = ""
        self.exp_file_str = ""
        self.dims = []
        self.MPFU = MPFU
        self.AWG = AWG
        self.pcyc_str = ""
        
        # Build hash of initial pulses
        for id, pulse in enumerate(sequence.pulses):
            if type(pulse) == autoEPR.Detection:
                self.var_hash[id] = "pg"
                self._addDef(f"pg = {int(pulse.tp.value)}")
            elif type(pulse) == autoEPR.Delay:
                str = self._new_delay(id)
                self._addDef(f"{str} = {int(pulse.tp.value)}")
            else:
                str = self._new_pulse(id)
                self._addDef(f"{str} = {int(pulse.tp.value)}")

        # Build table of parvars
        unique_parvars = np.unique(sequence.progTable["axID"])
        self.parvars = []
        for i in unique_parvars:
            parvar = PSparvar(sequence, i)
            self.parvars.append(parvar)
        
        self._addDef(f"h = {sequence.shots.value}")
        self._addDef(f"n = {sequence.averages.value}")
        self._addDef(f"srt = {sequence.reptime.value:.0f} * srtu")

        # TODO: Add checker that the only the first two parvars are bruker 
        # compatible.

        # Incremented variables need both a placeholder variable and a stepper.
        step_hash = {}
        place_hash = {}
        for i, parvar in enumerate(self.parvars):
            for j, event in enumerate(parvar.events): 
                str = self._new_delay(f"{event}_step")
                step_hash[event] = str
                self._addDef(f"{str} = {parvar.parvar['step'][j]}")
                str = self._new_delay(f"{event}_place")
                place_hash[event] = str
            dim_step = self._new_delay(f"dim{i}_step")
            self._addDef(f"{dim_step} = {parvar.ax_step}")
            self._addDef(f"d{20+i} = {parvar.init}")
            self.dims.append(parvar.dim)

        # Build experiment file

        if self.AWG:
            self.pcyc_hash = {}
            id = -1
            for i,pulse in enumerate(sequence.pulses):
                if (type(pulse) is autoEPR.Delay):
                    continue
                elif (type(pulse) is autoEPR.Detection):
                    continue
                id += 1
                awg_id = self._addAWGPulse(pulse, id)
                self.pcyc_hash[i] = awg_id

            pcyc = PSPhaseCycle(self.sequence, MPFU=None, OnlyDet=True)
            self.pcyc_str = pcyc.__str__() + self.pcyc_str
                
        else:
            self._addPhaseCycle()

        for i in place_hash.keys():
            self._addExp(f"{place_hash[i]} = {self.var_hash[i]}")
        if self.parvars[0].Bsweep:
            self._addExp(f"bsweep x=1 to sx")
        else:
            self._addExp(f"sweep x=1 to sx")
    
        self._addExp(f"shot i=1 to h")
        self._addExp(f"d9")  # Just a stupid uncessary variable

        # # Read pulse sequence with phase cycle
        for id, pulse in enumerate(sequence.pulses):
            if id in place_hash:
                pulse_str = place_hash[id]
            else:
                pulse_str = self.var_hash[id]
            if id in self.pcyc_hash:
                self._addExp(f"{pulse_str} [{self.pcyc_hash[id]}]")
            elif type(pulse) == autoEPR.Detection:
                self._addExp(f"d0\nacq [sg1]")
            else:
                self._addExp(f"{pulse_str}")
        # if type(pulse) == autoEPR.Delay:
        #     self._addExp(f"{}")
        # else:
        #     self._addExp(f"{} [{}]")

        self._addExp(f"next i")  # End of shots loop
        for i in place_hash.keys():
            self._addExp(f"{place_hash[i]} = {place_hash[i]} + {step_hash[i]}")
        if not self.parvars[0].Bsweep:
            self._addExp(f"dx = dx + {dim_step}")
        self._addExp(f"next x")  # End of scan loop

        self._addScanLoop()
        self._cmpl_Exp()
        self._cmpl_def()
        pass

    def _new_delay(self, key):
        str = self.possible_delays.pop(0)
        self.var_hash[key] = str
        return str
    
    def _new_pulse(self, key):
        str = self.possible_pulses.pop(0)
        self.var_hash[key] = str
        return str

    def save(self, filepath):
        fdef_path = f"{filepath}.def"
        fexp_path = f"{filepath}.exp"
        with open(fdef_path, "w") as fdef:
            fdef.write(header.format("definition"))
            fdef.write(self.def_file_str)
        with open(fexp_path, "w") as fexp:
            fexp.write(header.format("experiment"))
            fexp.write(self._ExpDefs())
            fexp.write("\n"*3)
            fexp.write(self.pcyc.__str__())
            fexp.write("\n"*3)
            fexp.write(self.exp_file_str)

        pass

    def _addDef(self, str):
        self.def_file_str += str
        self.def_file_str += "\n"

    def _addExp(self, str):
        self.exp_file_str += str
        self.exp_file_str += "\n"

    def _ExpDefs(self):
        str = "begin defs \n"
        if len(self.dims) == 1:
            str += f"dim s[{self.dims[0]}] \n"
        elif len(self.dims) == 2:
            str += f"dim s[{self.dims[0]},{self.dims[1]}] \n"
        str += "end defs \n"
        return str
        
    def _addScanLoop(self):
        self.exp_file_str = "for k=1 to n\n" + self.exp_file_str
        self.exp_file_str += "next k \n"
        pass

    def _addPhaseCycle(self):
        self.pcyc = PSPhaseCycle(self.sequence, self.MPFU)
        self.pcyc_hash = self.pcyc.pcyc_hash
        self.pcyc_str = self.pcyc.__str__()
    
    def _addAWGPulse(self, pulse, id):
        awg_id = id
        if type(pulse) is autoEPR.RectPulse:
            shape = 1
            init_freq = pulse.freq.value
            final_freq = init_freq
        

        def add_AWG_line(elements, comment=None, indent=2, repeat=1):
            string = ""
            string += " "*indent
            for i in range(0,repeat):
                if isinstance(elements,list) or isinstance(elements,np.ndarray):
                    for ele in elements:
                        string += f"{ele:<4}  "
                else:
                    string += f"{elements:<4}  "
            
            if comment is not None:
                string += ";  " + comment
            
            string += "\n"
            return string

        string = ""
        
        phase = np.round(np.degrees(pulse.pcyc["Phases"]))
        amplitude = 0
        num_cycles = len(phase)
        
        string += f"begin awg{awg_id}\n"
        string += add_AWG_line(
            init_freq, repeat=num_cycles, comment="Initial Frequency [MHz]")
        string += add_AWG_line(
            final_freq, repeat=num_cycles, comment="Final Frequency [MHz]")
        string += add_AWG_line(
            phase, repeat=1, comment="Phase")
        string += add_AWG_line(
            amplitude, repeat=num_cycles, comment="Amplitude")
        string += add_AWG_line(
            shape, repeat=num_cycles, comment="Shape")
        string += f"end awg{awg_id}\n"
    
        self.pcyc_str += string

        return f"awg{id}"
        



    def _cmpl_Exp(self):
        self.exp_file_str = f"begin exp \"{'auto'}\" [INTG QUAD]\n" +\
             self.exp_file_str
        self.exp_file_str += "end exp\n"

    def _cmpl_def(self):
        self.def_file_str = f"begin defs \n\n" +\
             self.def_file_str
        self.def_file_str += "end defs\n"

    def __str__(self) -> str:
        output_str = "DEF File \n" + "#"*79 + "\n"
        output_str += self.def_file_str
        output_str += "#"*79 + "\n"
        output_str += "EXP File \n" + "#"*79 + "\n"
        output_str += self._ExpDefs() + "\n \n" + self.pcyc_str + "\n \n" + self.exp_file_str

        return output_str

# =============================================================================
# Functions
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
    if os.path.isabs(ps_file[0]):
        base_path = ""
    else:
        base_path = MODULE_DIR

    if len(ps_file) == 1:
        # Assuming that both the EXP file and the DEF file have the same 
        # name bar-extention
        exp_file = base_path + ps_file[0] + ".exp"
        def_file = base_path + ps_file[0] + ".def"

    elif len(ps_file) == 2:
        
        # EXP and DEF file have seperate name
        exp_file = base_path + ps_file[0] + ".exp"
        def_file = base_path + ps_file[1] + ".def"

    else:
        raise ValueError(
            "ps_file must be of form ['EXP file'] or ['EXP file','DEF file']")

    # Identifying a dimension change in settings
    r = re.compile("dim([0-9]*)")
    match_list: list = list(filter(
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


def change_dimensions(path, dim: int, new_length):    
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


