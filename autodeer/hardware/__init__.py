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
    from .Bruker_tools import run_general,change_dimensions, PulseSpel
    from .XeprAPI_link import *
    from .Bruker_MPFU import *
    from .Bruker_AWG import *


