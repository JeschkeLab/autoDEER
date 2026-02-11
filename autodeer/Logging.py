import logging
import logging.handlers as handlers
import os
import PyQt6.QtCore as QtCore


class DictFormater(logging.Formatter):

    def format(self, record):
        
        entry = {
            'time': record.asctime,
            'name': record.name,
            'level': record.levelname,
            'message': record.msg
        }

        return entry
    

class QtLogHandler(QtCore.QObject, logging.Handler):

    signal = QtCore.pyqtSignal(dict)
    def __init__(self, level=logging.NOTSET):
        super().__init__(level=level)
        

    def emit(self, record):
        entry = self.format(record)
        self.signal.emit(entry)


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

    autoDEER_log = logging.getLogger('autoDEER')
    interface_log = logging.getLogger('interface')
    logHandler_core = handlers.TimedRotatingFileHandler(
        os.path.join(folder,'autoDEER.log'), when='D', backupCount=4)
    # logHandler_core = handlers.FileHandler(os.path.join(folder,'autoDEER.log'))
    logHandler_core.setFormatter(formatter)
    autoDEER_log.setLevel(logging.INFO)
    autoDEER_log.addHandler(logHandler_core)

    QTHandler = QtLogHandler()
    QTHandler.setFormatter(DictFormater())
    autoDEER_log.addHandler(QTHandler)

    logHandler_hardware = handlers.TimedRotatingFileHandler(
        os.path.join(folder,'interface.log'), when='D', backupCount=4)
    # logHandler_hardware = handlers.FileHandler(os.path.join(folder,'interface.log'))
    logHandler_hardware.setFormatter(formatter)
    interface_log.setLevel(logging.INFO)
    interface_log.addHandler(logHandler_hardware)
    interface_log.addHandler(QTHandler)


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
    logging.getLogger('autoDEER').info(f'Level set to {core}')
    logging.getLogger('interface').setLevel(interface)
    logging.getLogger('interface').info(f'Level set to {interface}')

