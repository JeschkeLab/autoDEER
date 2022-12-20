import autoDeer.hardware.openepr as autoEPR
import numpy as np
from autoDeer import __version__

header = "; " + "-" * 79 + "\n"
header += "; " + \
     f"Auto-generated PulseSpel {{}} file by autoDEER {__version__}\n"
header += "; " + "-" * 79 + "\n"


class PSPhaseCycle:

    def __init__(self, sequence) -> None:
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

        #  Detection event
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
        self.BPhaseCycles = BPhaseCycles

        self.__str__()

        pass

    def __str__(self) -> str:
        full_str = "begin lists \"auto\" \n"
        pcyc_hash = {}
        for i, cycle_dict in enumerate(self.BPhaseCycles):
            cycle = cycle_dict["Cycle"]
            if type(cycle_dict["Pulse"]) != str:
                pcyc_hash[cycle_dict["Pulse"]] = i+1
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


class PulseSpel:

    def __init__(self, sequence) -> None:
        self.possible_delays = possible_delays.copy()
        self.possible_pulses = possible_pulses.copy()
        self.sequence = sequence
        self.var_hash = {}
        self.def_file_str = ""
        self.exp_file_str = ""
        self.dims = []
        
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
        self._addPhaseCycle()

        for i in place_hash.keys():
            self._addExp(f"{place_hash[i]} = {self.var_hash[i]}")
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
                self._addExp(f"{pulse_str} [ph{self.pcyc_hash[id]}]")
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
        self._addExp(f"dx = dx + {dim_step}")
        self._addExp(f"next x")  # End of scan loop

        self._addScanLoop()
        self._cmpl_Exp()
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
        self.pcyc = PSPhaseCycle(self.sequence)
        self.pcyc_hash = self.pcyc.pcyc_hash

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
        output_str += self._ExpDefs() + "\n \n" + self.pcyc.__str__() + "\n \n" + self.exp_file_str

        return output_str
