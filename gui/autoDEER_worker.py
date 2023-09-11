import PyQt6.QtCore as QtCore
import sys, traceback
import autodeer as ad
import time
import numpy as np
from threadpoolctl import threadpool_limits


class autoDEERSignals(QtCore.QObject):
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
    progress = QtCore.pyqtSignal(int)
    status = QtCore.pyqtSignal(str)
    fsweep_result  = QtCore.pyqtSignal(object)
    respro_result = QtCore.pyqtSignal(object)
    optimise_pulses = QtCore.pyqtSignal(object)
    relax_result = QtCore.pyqtSignal(object)
    quickdeer_result = QtCore.pyqtSignal(object)
    longdeer_result = QtCore.pyqtSignal(object)
    

class autoDEERWorker(QtCore.QRunnable):
    '''
    Worker thread

    Inherits from QRunnable to handler worker thread setup, signals and wrap-up.

    :param callback: The function callback to run on this worker thread. Supplied args and
                     kwargs will be passed through to the runner.
    :type callback: function
    :param args: Arguments to pass to the callback function
    :param kwargs: Keywords to pass to the callback function

    '''

    def __init__(self, interface, wait:QtCore.QWaitCondition, mutex:QtCore.QMutex, results:dict, *args, **kwargs):
        super(autoDEERWorker,self).__init__()

        # Store constructor arguments (re-used for processing)
        self.interface = interface
        self.args = args
        self.kwargs = kwargs
        self.signals = autoDEERSignals()
        self.updaterate = 30
        self.wait = wait
        self.mutex = mutex
        self.results = results
        self.sample = self.kwargs['sample']


        # # Add the callback to our kwargs
        # self.kwargs['progress_callback'] = self.signals.progress
    def pause_and_wait(self):
        self.mutex.lock()
        self.wait.wait(self.mutex)
        self.mutex.unlock()


    def run_fsweep(self):
        '''
        Initialise the runner function with passed args, kwargs.
        '''

        LO = 33.98
        gyro_N = 0.0028087
        reptime = 3e3
        p90, p180 = self.interface.tune_rectpulse(tp=12, LO=LO, B=LO/gyro_N, reptime = reptime)

        fsweep = ad.FieldSweepSequence(
            B=LO/gyro_N, LO=LO,reptime=reptime,averages=1,shots=40,
            Bwidth = 500, 
            pi2_pulse=p90, pi_pulse=p180,
        )

        self.interface.launch(fsweep,savename=f"{self.sample}_fieldsweep")
        self.signals.status.emit('Field-sweep running')
        while self.interface.isrunning():
            time.sleep(self.updaterate)
        self.signals.status.emit('Field-sweep complete')
        self.signals.fsweep_result.emit(self.interface.acquire_dataset())

    def run_respro(self):
        '''
        Initialise the runner function for resonator profile.
        '''

        LO = 33.98
        gyro_N = 0.0028087
        reptime = 3e3
        p90, p180 = self.interface.tune_rectpulse(tp=12, LO=LO, B=LO/gyro_N, reptime = reptime)
        self.gyro_exp = self.results['fieldsweep'].gyro
        RPseq = ad.ResonatorProfileSequence(
            B=LO/self.gyro_exp, LO=LO,reptime=reptime,averages=1,shots=120,
            pi2_pulse=p90, pi_pulse=p180,
        )

        self.interface.launch(RPseq,savename=f"{self.sample}_resonator_profile",IFgain=2)
        self.signals.status.emit('Resonator Profile running')
        while self.interface.isrunning():
            time.sleep(self.updaterate)
        self.signals.status.emit('Resonator Profile complete')
        self.signals.respro_result.emit(self.interface.acquire_dataset())

    def run_relax(self):
        '''
        Initialise the runner function for relaxation. 
        '''

        LO = 33.98
        reptime = 3e3
        relax = ad.DEERSequence(
            B=LO/self.gyro_exp, LO=LO,reptime=reptime,averages=1,shots=150,
            tau1=0.5, tau2=0.5, tau3=0.2, dt=15,
            exc_pulse=self.pulses['exc_pulse'], ref_pulse=self.pulses['ref_pulse'],
            pump_pulse=self.pulses['pump_pulse'], det_event=self.pulses['det_event']
            )
        relax.five_pulse(relaxation=True, re_step=150)
        relax.select_pcyc("16step_5p")
        relax._estimate_time();
        relax.pulses[1].scale.value = 0
        relax.pulses[3].scale.value = 0
        self.interface.launch(relax,savename=f"{self.sample}_relaxation",IFgain=2)
        self.signals.status.emit('Relaxation experiment running')
        while self.interface.isrunning():
            time.sleep(self.updaterate)
        self.signals.status.emit('Relaxation experiment complete')
        self.signals.relax_result.emit(self.interface.acquire_dataset())

    def run_quick_deer(self):
        LO = 33.98
        reptime = 3e3
        deer = ad.DEERSequence(
            B=LO/self.gyro_exp, LO=LO,reptime=reptime,averages=1000,shots=150,
            tau1=2.0, tau2=2.0, tau3=0.3, dt=15,
            exc_pulse=self.pulses['exc_pulse'], ref_pulse=self.pulses['ref_pulse'],
            pump_pulse=self.pulses['pump_pulse'], det_event=self.pulses['det_event']
            )
        deer.five_pulse()
        deer.select_pcyc("16step_5p")
        deer._estimate_time();
        self.interface.launch(deer,savename=f"{self.sample}_quickdeer",IFgain=2)
        self.signals.status.emit('DEER experiment running')
        with threadpool_limits(limits=1, user_api='blas'):
            self.interface.terminate_at(ad.criteria.DEERCriteria(mode="speed",verbosity=2),verbosity=2,test_interval=0.5)
        self.signals.status.emit('DEER experiment complete')
        self.signals.quickdeer_result.emit(self.interface.acquire_dataset())

    def run_long_deer(self):
        LO = 33.98
        reptime = 3e3
        rec_tau = self.results['quickdeer'].rec_tau
        max_tau = self.results['relax'].max_tau
        tau = np.min([rec_tau,max_tau])
        deer = ad.DEERSequence(
            B=LO/self.gyro_exp, LO=LO,reptime=reptime,averages=1000,shots=150,
            tau1=tau, tau2=tau, tau3=0.3, dt=15,
            exc_pulse=self.pulses['exc_pulse'], ref_pulse=self.pulses['ref_pulse'],
            pump_pulse=self.pulses['pump_pulse'], det_event=self.pulses['det_event']
            )
        deer.five_pulse()
        deer.select_pcyc("16step_5p")
        deer._estimate_time();
        self.interface.launch(deer,savename=f"{self.sample}_deer_tau_{tau}us",IFgain=2)
        self.signals.status.emit('DEER experiment running')
        with threadpool_limits(limits=1, user_api='blas'):
            self.interface.terminate_at(ad.DEERCriteria(mode="high"),verbosity=2,test_interval=0.5)
        self.signals.status.emit('DEER experiment complete')
        self.signals.quickdeer_result.emit(self.interface.acquire_dataset())

    @QtCore.pyqtSlot()    
    def run(self):

        LO = 33.98
        gyro_N = 0.0028087
        reptime = 3e3
        p90, p180 = self.interface.tune_rectpulse(tp=12, LO=LO, B=LO/gyro_N, reptime = reptime)

        self.run_fsweep()
        
        self.pause_and_wait()

        self.run_respro()

        self.pause_and_wait()

        exc_pulse = ad.RectPulse(  
                tp=16, freq=0, flipangle=np.pi/2, scale=0
            )
        ref_pulse = ad.RectPulse(  
                        tp=16, freq=0, flipangle=np.pi, scale=0
                    )
        pump_pulse = ad.HSPulse(  
                        tp=120, init_freq=-0.25, final_freq=-0.03, flipangle=np.pi, scale=0,
                        order1=6, order2=1, beta=10
                    )
        det_event = ad.Detection(tp=512, freq=0)

        self.pulses = {'exc_pulse':exc_pulse, 'ref_pulse':ref_pulse, 'pump_pulse':pump_pulse, 'det_event':det_event}

        # Optimise the pulses
        self.signals.status.emit('Optimising pulses')
        self.signals.optimise_pulses.emit(self.pulses)
        self.pause_and_wait()

        self.signals.status.emit('Tuning pulses')
        # Tune the pulses
        exc_pulse = self.interface.tune_pulse(exc_pulse, mode="amp_nut", B=LO/gyro_N,LO=LO,reptime=reptime,shots=100)
        ref_pulse = self.interface.tune_pulse(ref_pulse, mode="amp_nut", B=LO/gyro_N,LO=LO,reptime=reptime,shots=100)
        pump_pulse = self.interface.tune_pulse(pump_pulse, mode="amp_nut", B=LO/gyro_N,LO=LO,reptime=reptime,shots=100)


        # Run a relaxation experiment
        self.signals.status.emit('Running relaxation Experiment')
        self.run_relax()
        self.signals.status.emit('Analysing relaxation Experiment')

        self.pause_and_wait()

        self.signals.status.emit('Running QuickDEER')
        self.run_quick_deer()
        self.signals.status.emit('Analysing QuickDEER')
        # self.pause_and_wait()

        # self.signals.status.emit('Running LongDEER')
        # self.run_longdeer()
        # self.pause_and_wait()

        self.signals.finished.emit()

    def new_data(self, data):
        self.results = data


