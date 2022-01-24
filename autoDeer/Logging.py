import logging
import logging.handlers as handlers
import time

def setup_logs():
    # Here we set up the logging format. This will be the same for every logger
    formatter = logging.Formatter('%(asctime)s [%(name)s] - %(levelname)s: %(message)s')

    # A new log is created every week, after 4 weeks the oldest is deleated.
    # This might be changed to daily at a later date, we will see how it goes
    logHandler_core = handlers.TimedRotatingFileHandler('core.log', when='W0', backupCount=4)
    logHandler_core.setLevel(logging.INFO)
    logHandler_core.setFormatter(formatter)
    logging.getLogger('core').addHandler(logHandler_core)


    logHandler_hardware = handlers.TimedRotatingFileHandler('hardware.log', when='W0', backupCount=4)
    logHandler_hardware.setLevel(logging.INFO)
    logHandler_hardware.setFormatter(formatter)
    logging.getLogger('hardware').addHandler(logHandler_hardware)


