
from PyQt6.QtWidgets import QApplication, QWidget,QLabel,QDoubleSpinBox,QGridLayout,QAbstractSpinBox
from PyQt6 import uic
import PyQt6.QtCore as QtCore
import PyQt6.QtGui as QtGui
from pathlib import Path
from threadpoolctl import threadpool_limits
import os

from matplotlib.backends.backend_qtagg import FigureCanvas, NavigationToolbar2QT
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
from autodeer import DEERanalysis, DEERanalysis_plot
import deerlab as dl

from autodeer.gui.tools import * 
from functools import partial
from scipy.stats import norm

QtCore.QDir.addSearchPath('icons', 'gui/resources')
package_directory = os.path.dirname(os.path.abspath(__file__))


def get_Vexp(dataset, tmin=0):
    # t = dataset.axes[0]
    t = dataset['t'].data
    Vexp = dataset.data

    # normalise v
    Vexp /= Vexp.max()

    if t.max()>200:
        # guess this is in us not ns
        t /= 1e3

    # t -= t.min()
    t += tmin

    return t, Vexp


def deeranalysis_process(dataset, settings, cores):

    with threadpool_limits(limits=cores, user_api='blas'):
        return DEERanalysis(dataset, **settings, verbosity=2)
    


class DEERplot(QWidget):
    def __init__(self,parent=None):
        
        super().__init__(parent)
 
        # loading the ui fsile with uic module
        uic.loadUi(f"{package_directory}/deerpanel.ui", self)

        self.threadpool = QtCore.QThreadPool()
        self.current_results = {}
        self.current_data = {}
        self.create_figure()
        self.setup_inputs()
        self.DL_params={}

        self.current_folder = ''
        self.cores = 1
        self.Last_updated.setText(f"Last updated: never")
        self.num_scans.setText(f"# of scans: 0")


    def setup_inputs(self):
        self.ExperimentcomboBox.addItems(['5pDEER', '4pDEER'])
        self.DistancecomboBox.addItems(['auto','Parametric', 'Gaussian', 'Bi-Gaussian', 'Tri-Gaussian','Rice','Bi-Rice','Tri-Rice'])
        self.BackgroundcomboBox.addItems(['Hom3D','Str-exp'])
        self.PathwayslineEdit.setText('1,2,3,4,5')
        self.CompactnessradioButton.setChecked(True)
        self.PulseLengthdoubleSpinBox.setValue(16)

    def update_inputs_from_dataset(self):
        dataset= self.current_data['quickdeer']
        if hasattr(dataset, 'sequence'):
            seq = getattr(dataset, 'sequence')
            self.ExperimentcomboBox.setCurrentText(seq.name)
            if seq.name =='4pDEER':
                self.Tau1doubleSpinBox.setValue(param_in_us(seq.tau1))
                self.Tau2doubleSpinBox.setValue(param_in_us(seq.tau2))
                self.Tau3doubleSpinBox.setDisabled(1)
                self.PathwayslineEdit.setText('1,2,3')
            if seq.name =='5pDEER':
                self.Tau3doubleSpinBox.setDisabled(0)

                self.Tau1doubleSpinBox.setValue(param_in_us(seq.tau1))
                self.Tau2doubleSpinBox.setValue(param_in_us(seq.tau2))
                self.Tau3doubleSpinBox.setValue(param_in_us(seq.tau3))
                self.PathwayslineEdit.setText('1,2,3,4,5')

        elif 't' in dataset.coords:
            self.Tau1doubleSpinBox.setValue( dataset.tau1*1e-3)
            self.Tau2doubleSpinBox.setValue( dataset.tau2*1e-3)
            if (dataset.attrs['seq_name'] == '5pDEER'):
                self.ExperimentcomboBox.setCurrentText('5pDEER')
                self.Tau3doubleSpinBox.setDisabled(0)
                self.Tau3doubleSpinBox.setValue( dataset.tau3*1e-3)
                self.PathwayslineEdit.setText('1,2,3,4,5')
            elif (dataset.attrs['seq_name'] == '4pDEER'):
                self.ExperimentcomboBox.setCurrentText('4pDEER')
                self.Tau3doubleSpinBox.setDisabled(1)
                self.PathwayslineEdit.setText('1,2,3')
            
        


    def clearLayout(self, layout):
            if layout is not None:
                while layout.count():
                    item = layout.takeAt(0)
                    widget = item.widget()
                    if widget is not None:
                        widget.deleteLater()
                    else:
                        self.clearLayout(item.layout())

    def update_analysis_table(self):
        results = self.fitresult
        headers= ['Parameter', 'Value', '95% CI', 'Unit']
        rows=[]

        # fitparams = {key : fitvalue if len(fitvalue)>1 else fitvalue[0] for key, fitvalue in zip(results.paramlist,[results.param[idx] for idx in results._param_idx])}

        for param in results.paramlist:
            if param == 'P':
                continue
            rows.append({"Parameter":param, "Value":getattr(results,param), "95% CI":getattr(results,f"{param}Uncert").ci(95), "Unit":""})
        fill_table(self.Analysis_table, headers, rows, rowcount = len(rows))

    def update_fit_result(self):

        results = self.fitresult

        self.MNRDoubleSpinBox.setValue(results.MNR)
        self.Chi2DoubleSpinBox.setValue(results.stats['chi2red'])
        self.conc_value.setValue(results.conc)
        self.conc_uncert.setText(getCIstring(results.concUncert))

        try:
            self.regparamDoubleSpinBox.setValue(results.regparam)
            self.regparamDoubleSpinBox.setDisabled(0)

        except AttributeError:
            self.regparamDoubleSpinBox.setDisabled(1)

        self.clearLayout(self.Pathways_Box)

        # Find pathways in the fit results
        pathways = []
        total_mod = 0
        total_mod_uncert_covmat = None
        total_mod_lb = 0
        total_mod_ub = 0
        if 'mod' in results.paramlist:
            pathways = [1]
            reftime = getattr(results,f"reftime")
            lam = getattr(results,f"mod")
            total_mod = lam
            total_mod_uncert_covmat = np.sqrt(getattr(results,f"modUncert").covmat)
            total_mod_lb = getattr(results,f"modUncert")._UQResult__lb
            total_mod_ub = getattr(results,f"modUncert")._UQResult__ub
            i=0
            self.Pathways_Box.addWidget(QLabel(f"reftime"),i,0,1,1)
            self.Pathways_Box.addWidget(QDoubleSpinBox(value=reftime, suffix=' us', readOnly=True, buttonSymbols=QAbstractSpinBox.ButtonSymbols(2)),i,1,1,1)
            self.Pathways_Box.addWidget(QLabel(getCIstring(getattr(results,f"reftimeUncert"))),i,2,1,1)
            self.Pathways_Box.addWidget(QLabel(f"mod"),i+1,0,1,1)
            self.Pathways_Box.addWidget(QDoubleSpinBox(value=lam, suffix='', decimals=3, readOnly=True, buttonSymbols=QAbstractSpinBox.ButtonSymbols(2)),i+1,1,1,1)
            self.Pathways_Box.addWidget(QLabel(getCIstring(getattr(results,f"modUncert"))),i+1,2,1,1)

        else:
            for param in results.paramlist:
                # Search for parameters that start with 'reftime' and record the number
                if param.startswith('reftime'):
                    pathways.append(int(param[7:]))
            
            # Add pathways to table
            for i in range(0,len(pathways)*2,2):
                pathway = pathways[i//2]

                # Creating a new lay out that is formed from a label and 2 double spin boxes and 2 more labels. in a 1x2x2 pattern
                reftime = getattr(results,f"reftime{pathway}")
                lam = getattr(results,f"lam{pathway}")
                if total_mod_uncert_covmat is not None:
                    total_mod_uncert_covmat += np.sqrt(getattr(results,f"lam{pathway}Uncert").covmat)
                else:
                    total_mod_uncert_covmat = np.sqrt(getattr(results,f"lam{pathway}Uncert").covmat)
                total_mod_lb += getattr(results,f"lam{pathway}Uncert")._UQResult__lb
                total_mod_ub += getattr(results,f"lam{pathway}Uncert")._UQResult__ub
                total_mod += lam
                self.Pathways_Box.addWidget(QLabel(f"reftime {pathway}"),i,0,1,1)
                self.Pathways_Box.addWidget(QDoubleSpinBox(value=reftime, suffix=' us', readOnly=True, buttonSymbols=QAbstractSpinBox.ButtonSymbols(2)),i,1,1,1)
                self.Pathways_Box.addWidget(QLabel(getCIstring(getattr(results,f"reftime{pathway}Uncert"))),i,2,1,1)
                self.Pathways_Box.addWidget(QLabel(f"lam {pathway}"),i+1,0,1,1)
                self.Pathways_Box.addWidget(QDoubleSpinBox(value=lam, suffix='', decimals=3, readOnly=True, buttonSymbols=QAbstractSpinBox.ButtonSymbols(2)),i+1,1,1,1)
                self.Pathways_Box.addWidget(QLabel(getCIstring(getattr(results,f"lam{pathway}Uncert"))),i+1,2,1,1)

        total_mod_uncert_covmat = np.sqrt(np.diag(total_mod_uncert_covmat))

        coverage = 95
        alpha = 1 - coverage/100
        p = 1 - alpha/2
        standardError = norm.ppf(p) * total_mod_uncert_covmat[0]
        minCI = np.max([total_mod - standardError,total_mod_lb[0]])
        maxCI = np.min([total_mod + standardError,total_mod_ub[0]])

        self.Total_mod_value.setValue(total_mod)
        self.Total_mod_uncert.setText(f"({minCI:.{2}f},{maxCI:.{2}f})")

    def create_figure(self):
        fig, axs = plt.subplot_mosaic([
            ['Primary_time', 'Primary_time', 'Primary_dist', 'Primary_dist']
            ], figsize=(12.5, 6.28))
        self.static_canvas = FigureCanvas(fig)
        Navbar = NavigationToolbar2QT(self.static_canvas, self)
        Navbar.setMaximumHeight(24)
        self.Plot_layout.addWidget(self.static_canvas)
        self.Plot_layout.addWidget(Navbar)

        self._static_ax = axs
        self.static_canvas.figure.tight_layout(pad=4)
        self.static_canvas.figure.subplots_adjust(bottom=0.2, hspace=0.4)
    
    
    def update_figure(self):
        if not hasattr(self,'fit_result'):
            # Only update the raw data
            if  hasattr(self,'current_data') and 'quickdeer' in self.current_data.keys():
                self._static_ax['Primary_time'].plot(*get_Vexp(self.current_data['quickdeer']))

        self.static_canvas.draw()

    def clear_all(self):
        self._static_ax['Primary_time'].cla()
        self._static_ax['Primary_dist'].cla()
        self.static_canvas.draw()
        self.Last_updated.setText(f"Last updated: never")
        self.num_scans.setText(f"# of scans: 0")

        self.MNRDoubleSpinBox.setValue(0)
        self.Chi2DoubleSpinBox.setValue(0)
        self.conc_value.setValue(0)
        self.conc_uncert.setText('(-,-)')
        self.regparamDoubleSpinBox.setValue(0)
        self.Total_mod_value.setValue(0)
        self.Total_mod_uncert.setText('(-,-)')
        # clear all enteries in the Pathways_Box
        while self.Pathways_Box.count() > 0:
            self.Pathways_Box.removeWidget(self.Pathways_Box.itemAt(0).widget())
        self.fitresult = None


        

    def process_deeranalysis(self,wait_condition=None, update_func=None):

        settings = {'ROI':True}
        settings['exp_type'] = self.ExperimentcomboBox.currentText()
        settings['tau1'] = self.Tau1doubleSpinBox.value()
        settings['tau2'] = self.Tau2doubleSpinBox.value()
        if settings['exp_type'] == '5pDEER':
            settings['tau3'] = self.Tau3doubleSpinBox.value()

        settings['pathways'] = str_to_list_type(self.PathwayslineEdit.text(), int)
        settings['compactness'] = self.CompactnessradioButton.isChecked()
        settings['pulselength'] = self.PulseLengthdoubleSpinBox.value()

        dataset= self.current_data['quickdeer']

        # if dataset.axes[0].max()>500: # Axis is in ns
        #     dataset.axes[0] /= 1e3
        
        # dataset.axes[0] -= dataset.axes[0].min()
        # dataset.axes[0] += self.tmindoubleSpinBox.value()

        if self.DistancecomboBox.currentText() == 'auto':
            settings['model'] = None
        elif self.DistancecomboBox.currentText() == 'Parametric':
            settings['model'] = None
        elif self.DistancecomboBox.currentText() == 'Gaussian':
            settings['model'] = dl.dd_gauss
        elif self.DistancecomboBox.currentText() == 'Bi-Gaussian':
            settings['model'] = dl.dd_gauss2
        elif self.DistancecomboBox.currentText() == 'Tri-Gaussian':
            settings['model'] = dl.dd_gauss3
        elif self.DistancecomboBox.currentText() == 'Rice':
            settings['model'] = dl.dd_rice
        elif self.DistancecomboBox.currentText() == 'Bi-Rice':
            settings['model'] = dl.dd_rice2
        elif self.DistancecomboBox.currentText() == 'Tri-Rice':
            settings['model'] = dl.dd_rice3
        else:
            settings['model'] = None


        worker = Worker(deeranalysis_process, dataset, settings, self.cores)
        worker.signals.result.connect(partial(self.refresh_deer, wait_condition=wait_condition,update_func=update_func))

        print('starting worker')
        self.threadpool.start(worker)

    def refresh_deer(self, fitresult, wait_condition=None,update_func=None):
        if isinstance(fitresult, tuple):
            self.fitresult = fitresult[0]
            self.taumax = fitresult[1]
        else:
            self.fitresult = fitresult
        if update_func is not None:
            update_func(self.fitresult)
        if isinstance(wait_condition , QtCore.QWaitCondition):
            wait_condition.wakeAll()
        

        self._static_ax['Primary_time'].cla()
        self._static_ax['Primary_dist'].cla()
        DEERanalysis_plot(self.fitresult, background=True, ROI=self.fitresult.ROI, axs= self._static_ax, fig=self.static_canvas.figure)
        self.static_canvas.draw()
        self.update_fit_result()

        if hasattr(self.fitresult,'dataset') and hasattr(self.fitresult.dataset,'time'):
            self.Last_updated.setText(f"Last updated: {self.fitresult.dataset.time}")
            self.num_scans.setText(f"# of scans: {self.fitresult.dataset.nAvgs}")
        else:
            self.Last_updated.setText(f"Last updated: never")
            self.num_scans.setText(f"# of scans: 0")


if __name__ == '__main__':
    app = QApplication([])
    window = DEERplot()
    window.show()
    app.exec()