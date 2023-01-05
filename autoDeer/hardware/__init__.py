from .keysight_awg import *
from .pulses import *
from .awg_experiments import deer_pulse_5p,sequence_nutation
from .openepr import *
try:
    import matlab.engine
except ModuleNotFoundError:
    pass
else:
    from .ETH_awg import *

try:
    import XeprAPI
except ModuleNotFoundError:
    pass
else:
    from .Bruker_tools import run_general,change_dimensions
    from .xepr_api_adv import *
    from .Bruker_MPFU import *


