from PyQt6.QtWidgets import QApplication, QMainWindow, QFileDialog,QTableWidgetItem, QMessageBox
from PyQt6 import uic
import PyQt6.QtCore as QtCore
import PyQt6.QtGui as QtGui
from pathlib import Path
import sys, traceback
import autodeer as ad
import numbers
import numpy as np

# =============================================================================
#
#   General functions and classes
#
# =============================================================================

def getCIstring(Uncert, precision=2):
    try:
        return f"({Uncert.ci(95)[0]:.{precision}f},{Uncert.ci(95)[1]:.{precision}f})"
    except ValueError:
        return "frozen"
    
def load_epr_file(Object, store_location):

        filename, _= QFileDialog.getOpenFileName(
            Object,"Select a File", Object.current_folder,"Data (*.DTA *.mat)")
        
        if filename:
                path = Path(filename)
                filename_edit = str(path)

        dataset = ad.eprload(filename_edit)
        Object.current_data[store_location] = dataset


def get_sequence_rows(Sequence, names:list):
     
    rows = []
    for param in names:
        attr = getattr(Sequence, param)
        if isinstance(attr, ad.Parameter):
            unit = attr.unit
        else:
            unit = None
        if isinstance(attr, ad.Parameter):
            rows.append({"Parameter":param, "Value":f"{attr.value:.2f}", "Unit":unit})
        else:
            rows.append({"Parameter":param, "Value":f"{attr:.2f}", "Unit":None})

def fill_table(table, headers, rows, rowcount=None):
     
    if rowcount is None:
        table.setRowCount(len(rows))
    else:
        table.setRowCount(rowcount)
    # Fill table from dictionary rows
    for i, row in enumerate(rows):
        for j, key in enumerate(headers):
            if key in row.keys():
                entry = row[key]
                if isinstance(entry, numbers.Number):
                    item = QTableWidgetItem(f"{entry:.3f}")
                elif isinstance(entry, np.ndarray):
                    item = QTableWidgetItem(np.array2string(entry, precision=3, separator=',',suppress_small=True,threshold=4))
                else:
                    item = QTableWidgetItem(str(row[key]))
                table.setItem(i, j, item)
            else:
                item = QTableWidgetItem('')
                table.setItem(i, j, item)


def read_table(table):
    rows = table.rowCount()
    columns = table.columnCount()
    headers = [str(table.horizontalHeaderItem(i).text()) for i in range(columns)]
    data = {}
    for row in range(rows):
        row_data = {}
        for col in range(columns):
            try:
                row_data[headers[col]] = table.item(row,col).data()
            except TypeError:
                row_data[headers[col]] = table.item(row,col).text()
            except AttributeError:
                if col ==0: 
                    break
                row_data[headers[col]] = None
        else:
            data[row_data[headers[0]]] = row_data

    return data

def list_str_to_type(list, type=int):
    new_list = []
    for item in list:
        new_list.append(type(item))

    return new_list
def str_to_list_type(string, type=int):
    list = string.split(',')
    new_list = []
    for item in list:
        new_list.append(type(item))

    return new_list

def pyqt_table_from_dict(table, dict):
    headers = list(dict[list(dict.keys())[0]].keys())
    rows = []
    for key in dict.keys():
        rows.append(dict[key])
    fill_table(table, headers, rows)

def param_in_us(param):
    if not isinstance(param, ad.Parameter):
        raise TypeError(f"Expected Parameter, got {type(param)}")
    if param.unit == 'us':
        return param.value
    elif param.unit == 'ns':
        return param.value/1000
    
def param_in_MHz(param):
    if not isinstance(param, ad.Parameter):
        raise TypeError(f"Expected Parameter, got {type(param)}")
    if param.unit == 'MHz':
        return param.value
    elif param.unit == 'GHz':
        return param.value*1000
    elif param.unit == 'kHz':
        return param.value/1000
    elif param.unit == 'Hz':
        return param.value/1000000
    
def test_SNR(Application, data):
    """Raises an error box if the SNR of the signal is less than 1.

    Parameters
    ----------
    data : _type_
        _description_
    """
    if data.epr.SNR < 1:
        QMessageBox.about(Application,'ERORR!', 'Signal to Noise ratio is less than 1. Please check the data and try again.')
        return False
    else:
        return True    
# =============================================================================
#
#   Multi-threading functions and classes
#
# =============================================================================

class WorkerSignals(QtCore.QObject):
    '''
    Defines the signals available from a running worker thread.

    Supported signals are:

    finished
        No data

    error
        tuple (exctype, value, traceback.format_exc() )

    result
        object data returned from processing, anything

    progress
        int indicating % progress

    '''
    finished = QtCore.pyqtSignal()
    error = QtCore.pyqtSignal(tuple)
    result = QtCore.pyqtSignal(object)

class Worker(QtCore.QRunnable):
    '''
    Worker thread

    Inherits from QRunnable to handler worker thread setup, signals and wrap-up.

    :param callback: The function callback to run on this worker thread. Supplied args and
                     kwargs will be passed through to the runner.
    :type callback: function
    :param args: Arguments to pass to the callback function
    :param kwargs: Keywords to pass to the callback function

    '''

    def __init__(self, fn, *args, **kwargs):
        super(Worker, self).__init__()

        # Store constructor arguments (re-used for processing)
        self.fn = fn
        self.args = args
        self.kwargs = kwargs
        self.signals = WorkerSignals()

        # # Add the callback to our kwargs
        # self.kwargs['progress_callback'] = self.signals.progress

    @QtCore.pyqtSlot()
    def run(self):
        '''
        Initialise the runner function with passed args, kwargs.
        '''

        # Retrieve args/kwargs here; and fire processing using them
        try:
            result = self.fn(*self.args, **self.kwargs)
        except:
            traceback.print_exc()
            exctype, value = sys.exc_info()[:2]
            self.signals.error.emit((exctype, value, traceback.format_exc()))
        else:
            self.signals.result.emit(result)  # Return the result of the processing
        finally:
            self.signals.finished.emit()  # Done

