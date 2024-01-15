from ._version import __version__
from .tools import eprload
from .utils import *
from .FieldSweep import FieldSweepAnalysis
from .ResPro import ResonatorProfileAnalysis, optimise_spectra_position
from .DEER_analysis import *
from .Relaxation import *
from .classes import *
from .sequences import *
from .pulses import *
from .criteria import *
from .dataset import *
from .reporter import Reporter

try:
    import PyQt6
except ImportError:
   pass 
else:
    from .gui.main import autoDEERUI
