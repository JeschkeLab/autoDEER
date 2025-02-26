from ._version import __version__
from .DEER_analysis import *
from .delay_optimise import *
from .sequences import *
from .criteria import *
from .reporter import Reporter, combo_figure
from .colors import primary_colors
from .pulse_optimise import build_default_pulses, create_pulses_rect, create_pulses_shape
from .delay_optimise import *
from .ref_2D_analysis import RefocusedEcho2DAnalysis

try:
    import PyQt6
except ImportError:
   pass 
else:
    from .gui.main import autoDEERUI
