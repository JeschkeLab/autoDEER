from autodeer.sequences import Sequence
from autodeer.pulses import RectPulse, Delay, Detection
from autodeer.utils import gcd, val_in_ns, val_in_us
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
            if type(pulse) == Delay:
                continue
            elif type(pulse) == Detection:
                continue
            dict = {}
            dict["Pulse"] = i
            dict["Cycle"] = []
            if pulse.pcyc["Channels"] == "ELDOR":
                dict["Cycle"].append("ELDOR")
            else:
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
            if sequence.pulses[i].pcyc["Channels"] == "ELDOR":
                dict["Cycle"].append("ELDOR")
            else:
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

    def __init__(self, sequence: Sequence, id) -> None:

        progTable = sequence.progTable
        progTable_n = len(progTable["EventID"])

        parvar = {}
        self.PulseSpel = True
        parvar["variables"] = []
        parvar["types"] = []
        self.events = []
        parvar["vec"] = []
        parvar["step"] = []

        for i in range(0, progTable_n):  # Loop over all progressions
            if progTable["axID"][i] == id:  # Ignore all progressions that aren't for this axis
                pulse_num = progTable["EventID"][i]
                var = progTable["Variable"][i]
                vec = progTable["axis"][i]
                parvar["variables"].append(var)
                if var == 'B':
                    continue
                
                if var == 't':
                    #add this to the assoicated delay
                    parvar["types"].append("d")
                    self.events.append(pulse_num)
                elif var == 'reptime':
                    parvar["types"].append('srt')
                    self.events.append('srt')
                elif not isinstance(sequence.pulses[pulse_num], Detection):
                    parvar["types"].append("p")
                    self.events.append(pulse_num)

                parvar["vec"].append(vec)

                if len(np.unique(np.diff(vec))) != 1:
                    self.PulseSpel = False
                    parvar["step"].append(0)
                else:
                    parvar["step"].append(np.unique(np.diff(vec))[0])
        
        if "B" in parvar["variables"]:
            self.Bsweep = True
        else:
            self.Bsweep = False

        if "reptime" in parvar["variables"]:
            self.repsweep = True
        else:
            self.repsweep = False
        
        if len(parvar['step']) == 0:
            self.ax_step = None
            self.init = None
            self.dim = None
            self.parvar = parvar
        else:
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

possible_vars = ['a','b','c','e']

possible_pulses = [
        "p0", "p1", "p2", "p3", "p4", "p5", "p6", "p7", "p8", "p9", "p10"]

def convert_progtable(progtable):
    """
    This function reformats the progtable to be compatible with the Bruker_tools. 
    This is done by converting the axis of each moving pulse to be relative to the previous moving pulse.
    """
    n_entries = len(progtable['EventID'])
    uaxes = np.unique(progtable['axID'])
    n_uaxes = len(uaxes)

    for ax in uaxes:
        base_axis = None
        total_axis = None
        for i in range(n_entries):
            if progtable['axID'][i] != ax:
                continue
            if progtable['Variable'][i] != 't':
                continue
            if base_axis is None:
                base_axis = progtable['axis'][i]
                total_axis = base_axis
            else:
                new_axis = progtable['axis'][i] - total_axis
                total_axis = progtable['axis'][i]
                progtable['axis'][i] = new_axis
    
    return progtable



def calc_rel_positions(sequence):
    """
    Calcuates the starting relative positions of all pulses in a sequence.
    """
    n_pulses = len(sequence.pulses)
    positions = np.zeros(n_pulses)
    for i,pulse in enumerate(sequence.pulses):
        positions[i] = pulse.t.value
    rel_positions = np.diff(positions,prepend=0)
    return rel_positions

# =============================================================================


class PulseSpel:

    def __init__(self, sequence, d0, MPFU=None, AWG=False) -> None:

        self._check_sequence(sequence)
        self.convert_progtable(sequence.progTable)
        self.possible_delays = possible_delays.copy()
        self.possible_pulses = possible_pulses.copy()
        self.possible_vars = possible_vars.copy()
        self.sequence = sequence
        self.var_hash = {}
        self.def_file_str = ""
        self.exp_file_str = ""
        self.dims = []
        self.MPFU = MPFU
        self.AWG = AWG
        self.pcyc_str = ""

        
        # Build hash of initial pulses
        # For every event and delay is created immediately prior to the event

        rel_positions = calc_rel_positions(sequence)
        for id, pulse in enumerate(sequence.pulses):
            str = self._new_delay(f"{id}d")
            self._addDef(f"{str} = {int(rel_positions[id])}") 
            
            if type(pulse) == Detection:
                self.var_hash[id] = "pg"
                self._addDef(f"pg = {int(pulse.tp.value)}")
            else:
                str = self._new_pulse(id)
                self._addDef(f"{str} = {int(pulse.tp.value)}")

        self._addDef(f"h = {sequence.shots.value}")
        self._addDef(f"n = {sequence.averages.value}")
        self._addDef(f"srt = {sequence.reptime.value:.0f} * srtu")
        self._addDef(f"d0 = {d0:.0f} ")


        # Build table of parvars
        unique_parvars = np.unique(sequence.progTable["axID"])
        self.parvars = []
        for i in unique_parvars:
            parvar = PSparvar(sequence, i)
            self.parvars.append(parvar)
        

        # TODO: Add checker that the only the first two parvars are bruker 
        # compatible.

        # Incremented variables need both a placeholder variable and a stepper.
        step_hash = {}
        place_hash = {}
        for i, parvar in enumerate(self.parvars):
            for j, event in enumerate(parvar.events): 
                if event == 'srt': #srt loop
                    str = self._new_var('srt')
                    self._addDef(f"{str} = {parvar.parvar['step'][j]} * srtu")
                    str = self._new_var('srt_step')
                    step_hash[event] = str
                    self._addDef(f"{str} = {parvar.parvar['step'][j]} * srtu")
                    place_hash[event] = 'srt'
                elif parvar.parvar['types'][j] == 'p': #pulse loop
                    str = self._new_delay(f"{event}_step")
                    step_hash[f"{event}p"] = str
                    self._addDef(f"{str} = {parvar.parvar['step'][j]}")
                    str = self._new_pulse(f"{event}_place")
                    place_hash[f"{event}p"] = str
                elif parvar.parvar['types'][j] == 'd': #delay loop
                    str = self._new_delay(f"{event}_step")
                    step_hash[f"{event}d"] = str
                    self._addDef(f"{str} = {parvar.parvar['step'][j]}")
                    str = self._new_delay(f"{event}_place")
                    place_hash[f"{event}d"] = str
            dim_step = self._new_delay(f"dim{i}_step")
            self._addDef(f"{dim_step} = {parvar.ax_step}")
            self._addDef(f"d{20+i} = {parvar.init}")
            self.dims.append(parvar.dim)

        # Build experiment file

        if self.AWG:
            self.pcyc_hash = {}
            id = -1
            for pulse_num,pulse in enumerate(sequence.pulses):
                if (type(pulse) is Delay):
                    continue
                elif (type(pulse) is Detection):
                    continue
                id += 1
                awg_id = self._addAWGPulse(sequence,pulse_num, id)
                self.pcyc_hash[pulse_num] = awg_id

            pcyc = PSPhaseCycle(self.sequence, MPFU=None, OnlyDet=True)
            self.pcyc_str = pcyc.__str__() + self.pcyc_str
                
        else:
            self._addPhaseCycle()

        for i in place_hash.keys():
            self._addExp(f"{place_hash[i]} = {self.var_hash[i]}")
        if self.parvars[0].Bsweep:
            self._addExp(f"bsweep x=1 to sx")
        elif self.parvars[0].repsweep:
            self._addExp(f"for x=1 to sx")
        else:
            self._addExp(f"sweep x=1 to sx")
    
        self._addExp(f"shot i=1 to h")
        self._addExp(f"d9")  # Just a stupid uncessary variable

        # # Read pulse sequence with phase cycle
        for id, pulse in enumerate(sequence.pulses):
            # if id in place_hash:
            #     pulse_str = place_hash[id]
            # else:
            #     pulse_str = self.var_hash[id]
            # Add delays
            if f"{id}d" in place_hash:
                delay_str = place_hash[f"{id}d"]
            else:
                delay_str = self.var_hash[f"{id}d"]
            self._addExp(f"{delay_str}")

            # Add pulses
            if f"{id}p" in place_hash:
                pulse_str = place_hash[f"{id}p"]
            else:
                pulse_str = self.var_hash[id]
            if id in self.pcyc_hash:
                self._addExp(f"{pulse_str} [{self.pcyc_hash[id]}]")
            elif type(pulse) == Detection:
                self._addExp(f"d0\nacq [sg1]")
            else:
                self._addExp(f"{pulse_str}")

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
    
    def _new_var(self, key):
        str = self.possible_vars.pop(0)
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
            fexp.write(self.pcyc_str)
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
    
    def _addAWGPulse(self, sequence, pulse_num, id):
        awg_id = id
        pulse=sequence.pulses[pulse_num]
        if type(pulse) is RectPulse:
            shape = 0
            init_freq = pulse.freq.value
            final_freq = init_freq
            if hasattr(pulse,"scale"):
                amplitude = pulse.scale.value
            else:
                amplitude = 0
        

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
        
        # phase = np.round(np.degrees(pulse.pcyc["Phases"]))
        phase = np.round(np.degrees(sequence.pcyc_cycles[:,sequence.pcyc_vars.index(pulse_num)]))

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
        

    def _check_sequence(self, sequence):

        for i, pulse in enumerate(sequence.pulses):
            if (type(pulse) is not Delay):
                continue
            elif (type(pulse) is not Detection):
                continue
            elif pulse.scale is None:
                raise RuntimeError(f"Missing scale in pulse {i}")

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

    if os.path.isfile(ps_file[0]):
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
        
        with open(exp_file, "r") as exp_file:
            exp_text = exp_file.read()

        with open(def_file, "r") as def_file:
            def_text = def_file.read()
    else:
        exp_text = ps_file[0]
        def_text = ps_file[1]

    # # Identifying a dimension change in settings
    # r = re.compile("dim([0-9]*)")
    # match_list: list = list(filter(
    #     lambda list: list is not None, [r.match(i) for i in variables.keys()]))
    # if len(match_list) >= 1:
    #     for i in range(0, len(match_list)):
    #         key = match_list[i][0]
    #         dim = int(r.findall(key)[0])
    #         new_length = int(variables[key])
    #         change_dimensions(exp_file, dim, new_length)
        


    verbMsgParam = api.cur_exp.getParam('*ftEPR.PlsSPELVerbMsg')
    plsSPELCmdParam = api.cur_exp.getParam('*ftEPR.PlsSPELCmd')
    api.XeprCmds.aqPgSelectBuf(1)

    api.cur_exp.getParam('*ftEpr.PlsSPELGlbTxt').value = def_text
    api.XeprCmds.aqPgShowDef()

    plsSPELCmdParam.value=3
    while not "The variable values are set up" in verbMsgParam.value:
        time.sleep(0.1)

    api.XeprCmds.aqPgSelectBuf(2)

    api.cur_exp.getParam('*ftEpr.PlsSPELPrgTxt').value = exp_text
    plsSPELCmdParam.value=7
    while not "Second pass ended" in verbMsgParam.value:
        time.sleep(0.1)


    if "ReplaceMode" in settings:
        api.set_ReplaceMode(settings["ReplaceMode"]) 
    else:
        api.set_ReplaceMode(False) 
    
    if "PhaseCycle" in settings:
        api.set_PhaseCycle(settings["PhaseCycle"]) 
    else:
        api.set_PhaseCycle(True) 

    if "Acquisition_mode" in settings:
        api.set_Acquisition_mode(settings["Acquisition_mode"])
    else:    
        api.set_Acquisition_mode(1)


    for var in variables:
        api.set_PulseSpel_var(var.lower(), variables[var])

    api.set_PulseSpel_experiment(exp[0])
    api.set_PulseSpel_phase_cycling(exp[1])

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

def _addAWGPulse(sequence, pulse_num, id, pcyc_str,amp_var=None):
    awg_id = id
    pulse=sequence.pulses[pulse_num]
    if type(pulse) is RectPulse:
        shape = 0
        init_freq = pulse.freq.value
        final_freq = init_freq
        if amp_var is not None:
            amplitude = amp_var
        else:
            if hasattr(pulse,"scale"):
                amplitude = pulse.scale.value
            else:
                amplitude = 0
        

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
    
    # phase = np.round(np.degrees(pulse.pcyc["Phases"]))
    phase = np.round(np.degrees(sequence.pcyc_cycles[:,sequence.pcyc_vars.index(pulse_num)]))

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

    pcyc_str += string

    return pcyc_str, f"awg{id}"

def get_arange(array):
    if np.ndim(array) > 1:
        raise ValueError("The array must be 1D")

    array = np.around(array,6)
    unique_diff = np.unique(np.diff(array))
    step = unique_diff[0]
    start = array[0]
    stop = array[-1]

    return start,stop+step, step


def build_unique_progtable(seq):
    progtable = seq.progTable
    unique_axs = np.unique(progtable["axID"])
    u_progtable = []
    d_shifts=[]
    for axID in unique_axs:
        axes = []
        variables = []
        for i,j in enumerate(progtable["axID"]):
            if axID != j:
                continue
            axes.append(progtable["axis"][i])
            variables.append((progtable["EventID"][i], progtable["Variable"][i]))
        
        if np.any([a[1]=='scale' for a in variables]):
            axes = [ax * 100 for ax in axes]

        # Find unique axis
        steps = []
        for axis in axes:
            start,_, step = get_arange(axis)
            steps.append(step)
        if np.any(steps):
            if np.all(np.mod(steps,1) == 0):
                common_step = gcd(steps)
            else:
                common_step = np.min(steps)
                # idx = np.argmin(steps)
                # start = get_arange(axes[idx])
            multipliers = np.array(steps)/common_step
        else:
            common_step = 0
            multipliers = np.ones(len(steps))

        start = axes[0][0]/multipliers[0]


        index = progtable["axID"].index(axID)
        common_axis={"start":start,"dim":np.shape(axes[0])[0],"step":common_step,"reduce":progtable["reduce"][index]}
        n_steps = len(steps)
        vars = [{"variable": variables[i],"multiplier":multipliers[i]} for i in range(n_steps)]
        n_pulses = len(seq.pulses)
        t_shift = np.zeros(n_pulses)
        for var,mul in zip(variables,multipliers):
            if var[1] != 't':
                continue
            t_shift[var[0]] = mul
        
        u_progtable.append({"axID":axID,"axis":common_axis,"variables":vars, "delay_shifts": np.diff(t_shift,prepend=0)})



    return u_progtable

# def _addAWGPulse(self, sequence, pulse_num, id):
#     awg_id = id
#     pulse=sequence.pulses[pulse_num]
#     if type(pulse) is RectPulse:
#         shape = 0
#         init_freq = pulse.freq.value
#         final_freq = init_freq
#         if hasattr(pulse,"scale"):
#             amplitude = pulse.scale.value
#         else:
#             amplitude = 0
    

#     def add_AWG_line(elements, comment=None, indent=2, repeat=1):
#         string = ""
#         string += " "*indent
#         for i in range(0,repeat):
#             if isinstance(elements,list) or isinstance(elements,np.ndarray):
#                 for ele in elements:
#                     string += f"{ele:<4}  "
#             else:
#                 string += f"{elements:<4}  "
        
#         if comment is not None:
#             string += ";  " + comment
        
#         string += "\n"
#         return string

#     string = ""
    
#     # phase = np.round(np.degrees(pulse.pcyc["Phases"]))
#     phase = np.round(np.degrees(sequence.pcyc_cycles[:,sequence.pcyc_vars.index(pulse_num)]))

#     num_cycles = len(phase)
    
#     string += f"begin awg{awg_id}\n"
#     string += add_AWG_line(
#         init_freq, repeat=num_cycles, comment="Initial Frequency [MHz]")
#     string += add_AWG_line(
#         final_freq, repeat=num_cycles, comment="Final Frequency [MHz]")
#     string += add_AWG_line(
#         phase, repeat=1, comment="Phase")
#     string += add_AWG_line(
#         amplitude, repeat=num_cycles, comment="Amplitude")
#     string += add_AWG_line(
#         shape, repeat=num_cycles, comment="Shape")
#     string += f"end awg{awg_id}\n"

#     self.pcyc_str += string

#     return f"awg{id}"

def check_variable(var:str, uprog):

    if isinstance(uprog,list):
        for up in uprog:
            if check_variable(var,up):
                return True
        return False
    else:
        for i in range(len(uprog['variables'])):
            if uprog['variables'][i]["variable"][1] == var:
                return True
        return False
    
def write_pulsespel_file(sequence,d0, AWG=False, MPFU=False):
    """Write the pulsespel file for a given sequence. 

    Parameters
    ----------
    sequence : Sequence
        The sequence class to be converted.
    AWG : bool, optional
        Is this a pulse spel file for an AWG spectrometer, by default False
    MPFU : list, optional
        A list of MPFU channels, by default False

    Returns
    -------
    str
        The string for the definition file
    str
        The string for the experiment file
    """

    uprogtable = build_unique_progtable(sequence)
    possible_delays = [
        "d1", "d2", "d3", "d4", "d5", "d6", "d7", "d8", "d10", "d11", "d12",
        "d13", "d14", "d15", "d16", "d17", "d18", "d19", "d20", "d21", "d22",
        "d23", "d24", "d25", "d26", "d27", "d28", "d29", "d30"]
    possible_pulses= ["p0", "p1","p2", "p3", "p4", "p5", "p6", "p7", "p8", "p9"]
    loop_iterators = ["i","j"]
    loop_dims = ["m","f"]
    letter_vars = ['b','c','e','f','g']
    possible_amps = ["aa1","aa2","aa3","aa4","aa5","aa6","aa7","aa8","aa9","aa10",
                     "ab1","ab2","ab3","ab4","ab5","ab6","ab7","ab8","ab9","ab10",
                     "ac1","ac2","ac3","ac4","ac5","ac6","ac7","ac8","ac9","ac10"]

    n_pulses = len(sequence.pulses)
    n_axes = len(uprogtable)

    static_delay_hash = {0:"d9"}
    static_pulse_hash = {}
    mv_delay_hash = {}
    mv_pulse_hash = {}
    axis_start_hash = {}
    axis_step_hash = {}
    axis_val_hash = {}
    mv_pulse_amp_hash = {}
    static_pulse_amp_hash = {}
    
    axs_ids = []


    head = "begin exp \"auto\" [INTG QUAD] \n"
    foot = "end exp"
    def_file = "begin defs \n"
    dims = "begin defs \n dim s["
    pcyc_str = ""

    # Add shot repetition time
    if not check_variable("reptime", uprogtable):
        def_file += f"srt = {val_in_us(sequence.reptime, False):.0f} * srtu\n"

    def_file += f"d0 = {d0:.0f}\n"
    prev_pulse = None
    for i,pulse in enumerate(sequence.pulses):
        static_delay_hash[i] = possible_delays.pop()
        if prev_pulse is None:
            def_file += f"{static_delay_hash[i]} = {pulse.t.value}\n"
        else:
            def_file += f"{static_delay_hash[i]} = {pulse.t.value - prev_pulse.t.value}\n"
        
        if type(pulse) is Detection:
            def_file += f"pg = {pulse.tp.value}\n"
        else:
            static_pulse_hash[i] = possible_pulses.pop()
            def_file += f"{static_pulse_hash[i]} = {pulse.tp.value}\n"

        if AWG and (not isinstance(pulse,Detection)):
            static_pulse_amp_hash[i] = possible_amps.pop()
            def_file += f"{static_pulse_amp_hash[i]} = {int(pulse.scale.value * 100):d}\n"
        prev_pulse = pulse

    for axis in uprogtable:
        ax_ID = axis["axID"]
        axs_ids.append(ax_ID)
        axis_start_hash[ax_ID] = possible_delays.pop()
        axis_step_hash[ax_ID] = possible_delays.pop()
        axis_val_hash[ax_ID] = possible_delays.pop()

        def_file += f"{axis_step_hash[ax_ID]} = {axis['axis']['step']}\n"
        if not axis['axis']['reduce']:
            dims += f"{axis['axis']['dim']},"

    axs_ids.sort(reverse=True)
    reduced_axes = np.array([uprogtable[i]["axis"]["reduce"] for i in range(n_axes)])
    n_reduced = np.sum(reduced_axes)

    if n_axes - n_reduced == 1:
        kept_axes_param = ["x"]
    elif n_axes - n_reduced == 2:
        kept_axes_param = ["x","y"]
    # elif n_axes - n_reduced == 0:
    #     raise RuntimeWarning("Transient experiments are not currently supported")
    # elif n_axes - n_reduced > 2:
    #     raise RuntimeWarning("A maximum of 2 non reduced axes are allowed.")



    # Setup averaging loop
    head += "for k=1 to m \ntotscans(m) \n"
    def_file += f"m = {sequence.averages.value}\n"

    foot = "scansdone(k) \nnext k \n"+foot

    delay_build = ""
    for ax in axs_ids:
        index = [uprogtable[i]['axID'] for i in range(len(uprogtable))].index(ax)
        uprog = uprogtable[index]

        if check_variable("B", uprog):

            reduce = False # You can not reduce a Bsweep
            loop_str = "bsweep {param}=1 to {stop} \n"
        else:
            reduce = uprogtable[index]['axis']["reduce"]
            if ax == 0 and not check_variable("reptime", uprog):
                loop_str = "sweep {param}=1 to {stop} \n"
            else:
                loop_str = "for {param}=1 to {stop} \n"
        
        if reduce:
            param = loop_iterators.pop()
        else:
            param=kept_axes_param.pop()

        foot = f"next {param} \n" + foot

        if param == 'x':
            stop = 'sx'
            foot = f"dx=dx+{axis_step_hash[ax]}\n"+ foot
        elif param == 'y':
            stop = 'sy'
            foot = f"dy=dy+{axis_step_hash[ax]}\n"+ foot
        else:
            stop = loop_dims.pop()
            def_file += f"{stop} = {uprog['axis']['dim']}\n"

        if check_variable("t", uprog):
            changing_pulses = uprog["delay_shifts"].nonzero()[0]
            for i in changing_pulses:
                if i not in mv_delay_hash:
                    mv_delay_hash[i] = possible_delays.pop()
                    # delay_build = f"{mv_delay_hash[i]} = {static_delay_hash[i]}\n" + delay_build 
                    head += f"{mv_delay_hash[i]} = {static_delay_hash[i]}\n"
                for k in range(int(np.abs(uprog["delay_shifts"][i]))):
                    if uprog["delay_shifts"][i] > 0: 
                        sign = "+"
                    else: 
                        sign = "-"
                    # foot = f"{axis_val_hash[ax]}={axis_val_hash[ax]}{sign}{axis_step_hash[ax]}\n"+ foot
                    foot = f"{mv_delay_hash[i]}={mv_delay_hash[i]}{sign}{axis_step_hash[ax]}\n"+ foot
                # delay_build += f"{mv_delay_hash[i]} = {mv_delay_hash[i]} + {axis_val_hash[ax]} \n"
        if check_variable("tp", uprog):
            for var in uprog["variables"]:
                if var["variable"][1] != "tp":
                    continue
                pulse_id = var["variable"][0]
                mult = var["multiplier"]
                mv_pulse_hash[pulse_id] = possible_pulses.pop()
                head += f"{mv_pulse_hash[pulse_id]} = {static_pulse_hash[pulse_id]}\n"
                if mult > 0:
                    sign = "+"
                else:
                    sign = "-"
                for k in range(int(np.abs(mult))):
                    foot = f"{mv_pulse_hash[pulse_id]}={mv_pulse_hash[pulse_id]}{sign}{axis_step_hash[ax]}\n"+ foot
        if check_variable("reptime", uprog):
            for var in uprog["variables"]:
                if var["variable"][1] != "reptime":
                    continue
                reptime0_var = letter_vars.pop()
                head += f"srt = {reptime0_var}\n"
                axis = val_in_us(sequence.reptime, True)
                def_file += f"{reptime0_var} = {axis[0]:.0f} * srtu\n"
                srt_step_var = letter_vars.pop()
                def_file += f"{srt_step_var} = {var['multiplier']*uprog['axis']['step']:.0f} * srtu\n"
                foot = f"srt=srt + {srt_step_var}\n"+ foot
        if check_variable("scale", uprog):
            for var in uprog["variables"]:
                if var["variable"][1] != "scale":
                    continue
                pulse_id = var["variable"][0]
                mult = var["multiplier"]
                mv_pulse_amp_hash[pulse_id] = possible_amps.pop()
                head += f"{mv_pulse_amp_hash[pulse_id]} = {static_pulse_amp_hash[pulse_id]}\n"
                
                if mult > 0:
                    sign = "+"
                else:
                    sign = "-"
                for k in range(int(np.abs(mult))):
                    foot = f"{mv_pulse_amp_hash[pulse_id]}={mv_pulse_amp_hash[pulse_id]}{sign}{axis_step_hash[ax]}\n"+ foot


        # head +=  f"{axis_val_hash[ax]} = 0\n"
        head += loop_str.format(param=param,stop=stop)



    # Add moving pulses
    # Add shot loop
    head += delay_build
    head += f"shot i=1 to h\n"
    def_file += f"h = {sequence.shots.value}\n"

    foot = "next i\n" + foot

    if AWG:
        pcyc_hash = {}
        id = -1
        for pulse_num,pulse in enumerate(sequence.pulses):
            if (type(pulse) is Detection):
                continue
            id += 1
            amp_var = mv_pulse_amp_hash[pulse_num] if pulse_num in mv_pulse_amp_hash else static_pulse_amp_hash[pulse_num]
            pcyc_str, awg_id = _addAWGPulse(sequence,pulse_num, id, pcyc_str,amp_var)
            pcyc_hash[pulse_num] = awg_id

        pcyc = PSPhaseCycle(sequence, MPFU=None, OnlyDet=True)
        pcyc_str = pcyc.__str__() + "\n"*1 + pcyc_str
    elif MPFU is not None:
        pcyc = PSPhaseCycle(sequence, MPFU)
        pcyc_hash = pcyc.pcyc_hash
        pcyc_str = pcyc.__str__()
    else:
        pcyc = PSPhaseCycle(sequence, MPFU)
        pcyc_hash = pcyc.pcyc_hash
        pcyc_str = pcyc.__str__()
    # Add Pulse Sequence
    pulse_str = "d9\n"
    for id, pulse in enumerate(sequence.pulses):
        if id in mv_delay_hash:
            pulse_str += f"{mv_delay_hash[id]}\n"
        else:
            pulse_str += f"{static_delay_hash[id]}\n"
        
        if type(pulse) == Detection:
            pulse_str += f"d0\nacq [sg1]\n"
        else:
            if id in mv_pulse_hash:
                pulse_str += f"{mv_pulse_hash[id]} [{pcyc_hash[id]}]\n"
            else:
                pulse_str += f"{static_pulse_hash[id]} [{pcyc_hash[id]}]\n"


    dims = dims[:-1]
    dims += "] \nend defs\n" 

    def_file += "end defs"
    exp_file = dims + "\n"*2 + pcyc_str + "\n"*2 +  head + pulse_str + foot

    return def_file, exp_file