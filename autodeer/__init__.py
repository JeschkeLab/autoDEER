from ._version import __version__
from .tools import eprload
from .utils import *
from .FieldSweep import FieldSweepAnalysis
from .ResPro import ResonatorProfileAnalysis, optimise_spectra_position
from .DEER_analysis import DEERanalysis,DEERanalysis_plot, DEERanalysis_plot_pub,IdentifyROI,optimise_pulses,plot_overlap,normalise_01, calc_DEER_settings
from .Relaxation import *
from .relaxation_autodeer import *
from .classes import Parameter, Interface
from .sequences import *
from .pulses import *
from .criteria import *
from .dataset import *
from .reporter import Reporter, combo_figure
from .colors import primary_colors

try:
    import PyQt6
except ImportError:
   pass 
else:
    from .gui.main import autoDEERUI
