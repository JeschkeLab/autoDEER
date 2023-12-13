import logging
import logging.handlers as handlers
import os

def setup_logs(folder: str):
    """
    General command to setup system wide logging.

    Two logs are setup, a core log and a harware log. The hardware logs records
    everything releated to the hardware interfaces whilst the core logs records
    all the analysis and related modules. 
    """
    # Here we set up the logging format. This will be the same for every logger
    formatter = logging.Formatter(
        '%(asctime)s [%(name)s] - %(levelname)s: %(message)s')

    # A new log is created every week, after 4 weeks the oldest is deleated.
    # This might be changed to daily at a later date, we will see how it goes
    logHandler_core = handlers.TimedRotatingFileHandler(
        os.path.join(folder,'autoDEER.log'), when='D', backupCount=4)
    logHandler_core.setLevel(logging.INFO)
    logHandler_core.setFormatter(formatter)
    logging.getLogger('autoDEER').addHandler(logHandler_core)

    logHandler_hardware = handlers.TimedRotatingFileHandler(
        os.path.join(folder,'interface.log'), when='D', backupCount=4)
    logHandler_hardware.setLevel(logging.INFO)
    logHandler_hardware.setFormatter(formatter)
    logging.getLogger('interface').addHandler(logHandler_hardware)

def change_log_level(core='INFO',interface='INFO'):
    """
    Change the log level of the core and hardware loggers.

    Parameters
    ----------
    core : str
        The log level for the core logger.
    interface : str
        The log level for the hardware logger.
    """
    logging.getLogger('autoDEER').setLevel(core)
    logging.getLogger('interface').setLevel(interface)
