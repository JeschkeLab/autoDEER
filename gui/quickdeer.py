
from PyQt6.QtWidgets import QApplication, QMainWindow, QFileDialog,QTableWidgetItem, QWidget, QVBoxLayout
from PyQt6 import uic
import PyQt6.QtCore as QtCore
import PyQt6.QtGui as QtGui
from pathlib import Path
import sys, traceback
from threadpoolctl import threadpool_limits

from matplotlib.backends.backend_qtagg import FigureCanvas, NavigationToolbar2QT
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import autodeer as ad
import deerlab as dl

import tools
from functools import partial

QtCore.QDir.addSearchPath('icons', '/Users/hugo/Documents/Work/autoDeer_dev/GUI/resources')


def get_Vexp(dataset, settings):
    t = dataset.axes[0]
    Vexp = dataset.data

    # normalise v
    Vexp /= Vexp.max()

    if t.max()>200:
        # guess this is in us not ns
        t /= 1e3

    if 'zerotime' in settings:
        t -= float(settings['zerotime']['Value'])

    return t, Vexp


def deeranalysis_process(dataset, settings):

    with threadpool_limits(limits=1, user_api='blas'):
        return ad.DEERanalysis(dataset, **settings, verbosity=2)
    


class DEERplot(QWidget):
    def __init__(self,parent=None):
        
        super().__init__(parent)
 
        # loading the ui fsile with uic module
        uic.loadUi("/Users/hugo/Documents/Work/autoDeer_dev/GUI/quickdeer.ui", self)

        self.threadpool = QtCore.QThreadPool()
        self.current_results = {}
        self.current_data = {}
        self.create_figure()
        self.toolbar()

        self.current_folder = ''


    def toolbar(self):
        upload_icon = QtGui.QIcon('icons:upload.png')

        def custom_load():
            tools.load_epr_file(self, 'quickdeer')
            self.update_exp_table()
            self.update_figure()
            

        self.Loadfile.setIcon(upload_icon)
        self.Loadfile.clicked.connect(custom_load)

        refresh_icon = QtGui.QIcon('icons:refresh.png')
        self.Refresh_analysis.setIcon(refresh_icon)
        self.Refresh_analysis.clicked.connect(self.process_deeranalysis)

    def update_exp_table(self):
        dataset= self.current_data['quickdeer']
        headers= ['Parameter', 'Value', 'Unit']
        rows=[]
        if hasattr(dataset, 'sequence'):
            seq = getattr(dataset, 'sequence')
            rows.append({"Parameter":"Type", "Value":seq.name, "Unit":""})
            if seq.name =='4pDEER':
                rows.append({"Parameter":"tau1", "Value":seq.tau1.value, "Unit":seq.tau1.unit})
                rows.append({"Parameter":"tau2", "Value":seq.tau2.value, "Unit":seq.tau2.unit})
            if seq.name =='5pDEER':
                rows.append({"Parameter":"tau1", "Value":seq.tau1.value, "Unit":seq.tau1.unit})
                rows.append({"Parameter":"tau2", "Value":seq.tau2.value, "Unit":seq.tau2.unit})
                rows.append({"Parameter":"tau3", "Value":seq.tau3.value, "Unit":seq.tau3.unit})
        else:
            rows.append({"Parameter":"Type", "Value":"5pDEER", "Unit":""})
            rows.append({"Parameter":"tau1", "Value":1.2, "Unit":'us'})
            rows.append({"Parameter":"tau2", "Value":1.2, "Unit":'us'})
            rows.append({"Parameter":"tau3", "Value":"0.3", "Unit":"us"})
            rows.append({"Parameter":"zerotime", "Value":"0.5", "Unit":"us"})

        rows.append({"Parameter":"Compactness", "Value":"True", "Unit":""})
        rows.append({"Parameter":"Pathways", "Value":"1,2,3,4,5", "Unit":""})
        tools.fill_table(self.Experiment_table, headers, rows, rowcount = 10)

    def update_analysis_table(self):
        results = self.fitresult
        headers= ['Parameter', 'Value', '95% CI', 'Unit']
        rows=[]

        # fitparams = {key : fitvalue if len(fitvalue)>1 else fitvalue[0] for key, fitvalue in zip(results.paramlist,[results.param[idx] for idx in results._param_idx])}

        for param in results.paramlist:
            if param == 'P':
                continue
            rows.append({"Parameter":param, "Value":getattr(results,param), "95% CI":getattr(results,f"{param}Uncert").ci(95), "Unit":""})
        tools.fill_table(self.Analysis_table, headers, rows, rowcount = len(rows))

    def get_exp_table(self):
        dict = tools.read_table(self.Experiment_table)
        return dict

    def create_figure(self):
        fig, axs = plt.subplot_mosaic([
            ['Primary_time', 'Primary_time', 'Primary_dist', 'Primary_dist']
            ], figsize=(12.5, 6.28))
        self.static_canvas = FigureCanvas(fig)
        self.Plot_layout.addWidget(NavigationToolbar2QT(self.static_canvas, self))
        self.Plot_layout.addWidget(self.static_canvas)
        # self._static_ax = self.static_canvas.figure.subplots(1,2)
        self._static_ax = axs
        self.static_canvas.figure.tight_layout(pad=4)
        self.static_canvas.figure.subplots_adjust(bottom=0.2, hspace=0.4)
    
    
    def update_figure(self):
        if not hasattr(self,'fit_result'):
            # Only update the raw data
            if  hasattr(self,'current_data') and 'quickdeer' in self.current_data.keys():
                self._static_ax['Primary_time'].plot(*get_Vexp(self.current_data['quickdeer'],self.get_exp_table()))

        self.static_canvas.draw()

        

    def process_deeranalysis(self,wait_condition=None):

        settings = {'ROI':True}
        table = self.get_exp_table()
        settings['exp_type'] = table['Type']['Value']
        settings['tau1'] = float(table['tau1']['Value'])
        settings['tau2'] = float(table['tau2']['Value'])
        if settings['exp_type'] == '5pDEER':
            settings['tau3'] = float(table['tau3']['Value'])

        settings['pathways'] = tools.str_to_list_type(table['Pathways']['Value'])
        settings['compactness'] = bool(table['Compactness']['Value'])
        # settings['compactness'] = False
        dataset= self.current_data['quickdeer']

        if dataset.axes[0].max()>500: # Axis is in ns
            dataset.axes[0] /= 1e3
        


        worker = tools.Worker(deeranalysis_process, dataset, settings)
        worker.signals.result.connect(partial(self.refresh_deer, wait_condition=wait_condition))

        print('starting worker')
        self.threadpool.start(worker)

    def refresh_deer(self, fitresult, wait_condition=None):
        if isinstance(fitresult, tuple):
            self.fitresult = fitresult[0]
            self.taumax = fitresult[1]
        else:
            self.fitresult = fitresult
        
        if isinstance(wait_condition , QtCore.QWaitCondition):
            wait_condition.wakeAll()

        self._static_ax['Primary_time'].cla()
        self._static_ax['Primary_dist'].cla()
        ad.DEERanalysis_plot(self.fitresult, background=False, ROI=self.fitresult.ROI, axs= self._static_ax, fig=self.static_canvas.figure)
        self.static_canvas.draw()
        self.update_analysis_table()


if __name__ == '__main__':
    app = QApplication([])
    window = DEERplot()
    window.show()
    app.exec()