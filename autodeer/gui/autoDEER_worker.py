import PyQt6.QtCore as QtCore
from autodeer import FieldSweepSequence, ResonatorProfileSequence, DEERSequence, RectPulse, ChirpPulse, HSPulse, Detection, DEERCriteria, SNRCriteria, ReptimeScan
import time
import numpy as np
from threadpoolctl import threadpool_limits


class autoDEERSignals(QtCore.QObject):
    '''
    Defines the signals available from a running worker thre

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
    quickdeer_update = QtCore.pyqtSignal(object)
    longdeer_update = QtCore.pyqtSignal(object)
    reptime_scan_result = QtCore.pyqtSignal(object)
    

class autoDEERWorker(QtCore.QRunnable):
    '''
    Worker thread

    Inherits from QRunnable to handler worker thread setup, signals and wrap-up.

    :param callback: The function callback to run on this worker thre Supplied args and
                     kwargs will be passed through to the runner.
    :type callback: function
    :param args: Arguments to pass to the callback function
    :param kwargs: Keywords to pass to the callback function

    '''

    def __init__(self, interface, wait:QtCore.QWaitCondition, mutex:QtCore.QMutex, results:dict, LO, gyro, user_inputs:dict = None, *args, **kwargs):
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
        self.samplename = user_inputs['sample']
        self.LO = LO
        self.gyro = gyro
        self.user_inputs = user_inputs

        if not 'DEER_update_func' in self.user_inputs:
            self.user_inputs['DEER_update_func'] = None
        
        if 'cores' in kwargs:
            self.cores = kwargs['cores']
        else:
            self.cores = 1


        # # Add the callback to our kwargs
        # self.kwargs['progress_callback'] = self.signals.progress
    def pause_and_wait(self):
        self.mutex.lock()
        self.wait.wait(self.mutex)
        self.mutex.unlock()


    def run_fsweep(self, reptime):
        '''
        Initialise the runner function with passed args, kwargs.
        '''
        LO = self.LO
        gyro_N = self.gyro
        p90, p180 = self.interface.tune_rectpulse(tp=12, LO=LO, B=LO/gyro_N, reptime = reptime,shots=100)

        fsweep = FieldSweepSequence(
            B=LO/gyro_N, LO=LO,reptime=reptime,averages=1,shots=40,
            Bwidth = 500, 
            pi2_pulse=p90, pi_pulse=p180,
        )

        self.interface.launch(fsweep,savename=f"{self.samplename}_fieldsweep")
        self.signals.status.emit('Field-sweep running')
        while self.interface.isrunning():
            time.sleep(self.updaterate)
        self.signals.status.emit('Field-sweep complete')
        self.signals.fsweep_result.emit(self.interface.acquire_dataset())

    def run_respro(self, reptime):
        '''
        Initialise the runner function for resonator profile.
        '''
        LO = self.LO
        gyro=self.gyro
        p90, p180 = self.interface.tune_rectpulse(tp=12, LO=LO, B=LO/gyro, reptime = reptime,shots=100)

        RPseq = ResonatorProfileSequence(
            B=LO/gyro, LO=LO,reptime=reptime,averages=5,shots=40,
            pi2_pulse=p90, pi_pulse=p180,
        )

        self.interface.launch(RPseq,savename=f"{self.samplename}_resonator_profile",IFgain=2)
        self.signals.status.emit('Resonator Profile running')

        self.interface.terminate_at(SNRCriteria(5))
        # while self.interface.isrunning():
        #     time.sleep(self.updaterate)
        self.signals.status.emit('Resonator Profile complete')
        self.signals.respro_result.emit(self.interface.acquire_dataset())

    def run_relax(self,reptime):
        '''
        Initialise the runner function for relaxation. 
        '''
        LO = self.LO
        gyro = self.gyro
        relax = DEERSequence(
            B=LO/gyro, LO=LO,reptime=reptime,averages=10,shots=20,
            tau1=0.5, tau2=0.5, tau3=0.2, dt=15,
            exc_pulse=self.pulses['exc_pulse'], ref_pulse=self.pulses['ref_pulse'],
            pump_pulse=self.pulses['pump_pulse'], det_event=self.pulses['det_event']
            )
        relax.five_pulse(relaxation=True, re_step=200)
        relax.select_pcyc("16step_5p")
        relax._estimate_time();
        relax.pulses[1].scale.value = 0
        relax.pulses[3].scale.value = 0
        self.interface.launch(relax,savename=f"{self.samplename}_relaxation",IFgain=2)
        self.signals.status.emit('Relaxation experiment running')
        self.interface.terminate_at(SNRCriteria(10))
        # while self.interface.isrunning():
        #     time.sleep(self.updaterate)
        self.signals.status.emit('Relaxation experiment complete')
        self.signals.relax_result.emit(self.interface.acquire_dataset())

    def run_quick_deer(self,reptime,tau=2.0, dt=15):
        LO = self.LO
        gyro = self.gyro
        deer = DEERSequence(
            B=LO/gyro, LO=LO,reptime=reptime,averages=50,shots=15,
            tau1=tau, tau2=tau, tau3=0.3, dt=dt,
            exc_pulse=self.pulses['exc_pulse'], ref_pulse=self.pulses['ref_pulse'],
            pump_pulse=self.pulses['pump_pulse'], det_event=self.pulses['det_event']
            )
        
        deer.five_pulse()
        deer.select_pcyc("16step_5p")
        deer._estimate_time();
        self.interface.launch(deer,savename=f"{self.samplename}_quickdeer",IFgain=2)
        self.signals.status.emit('DEER experiment running')
        with threadpool_limits(limits=self.cores, user_api='blas'):
            self.interface.terminate_at(DEERCriteria(mode="speed",verbosity=2,update_func=self.signals.quickdeer_update.emit),verbosity=2,test_interval=0.5)
        self.signals.status.emit('DEER experiment complete')
        self.signals.quickdeer_result.emit(self.interface.acquire_dataset())

    def run_long_deer(self, reptime,tau1,tau2,tau3=0, dt=15, deertype='5pDEER'):
        LO = self.LO
        
        
        deer = DEERSequence(
            B=LO/self.gyro, LO=LO,reptime=reptime,averages=500,shots=45,
            tau1=tau1, tau2=tau2, tau3=tau3, dt=dt,
            exc_pulse=self.pulses['exc_pulse'], ref_pulse=self.pulses['ref_pulse'],
            pump_pulse=self.pulses['pump_pulse'], det_event=self.pulses['det_event']
            )
        
        if deertype == '5pDEER':
            deer.five_pulse()
            deer.select_pcyc("16step_5p")
            savename_suffix = f"5p_{tau1:.3f}us_{tau2:.3f}us_{tau3:.3f}us"
        elif deertype == '4pDEER':
            deer.four_pulse()
            deer.select_pcyc("16step_4p")
            savename_suffix = f"4p_{tau1:.3f}us_{tau2:.3f}us"

        elif deertype == '3pDEER':
            deer.four_pulse()
            deer.select_pcyc("8step_3p")
            savename_suffix = f"3p_{tau1:.3f}u"

        deer._estimate_time();
        self.interface.launch(deer,savename=f"{self.samplename}_deer"+savename_suffix,IFgain=2)
        self.signals.status.emit('DEER experiment running')
        time.sleep(30) # Always wait for the experiment to properly start
        with threadpool_limits(limits=self.cores, user_api='blas'):
            self.interface.terminate_at(DEERCriteria(mode="high",verbosity=2,update_func=self.signals.longdeer_update.emit),verbosity=2,test_interval=0.5)
        self.signals.status.emit('DEER experiment complete')
        self.signals.longdeer_result.emit(self.interface.acquire_dataset())

    def run_reptime_opt(self,reptime_guess):
        LO = self.LO
        p90, p180 = self.interface.tune_rectpulse(tp=12, LO=LO, B=LO/self.gyro, reptime = reptime_guess,shots=100)

        scan = ReptimeScan(B=LO/self.gyro, LO=LO,reptime=reptime_guess, reptime_max=12e3, averages=1, shots=50,
                           pi2_pulse=p90, pi_pulse=p180)
        self.interface.launch(scan,savename=f"{self.samplename}_reptimescan",IFgain=2)
        # self.interface.terminate_at(SNRCriteria(15),verbosity=2,test_interval=0.5)
        self.signals.status.emit('Reptime scan complete')
        self.signals.reptime_scan_result.emit(self.interface.acquire_dataset())



    @QtCore.pyqtSlot()    
    def run(self):

        reptime = 3e3
        LO = self.LO
        gyro_N = self.gyro

        self.run_fsweep(reptime)
        
        self.pause_and_wait()

        self.run_reptime_opt(reptime)

        self.pause_and_wait()
        
        reptime = self.results['reptime'].optimal
        
        self.run_respro(reptime)
        
        
        self.pause_and_wait()

        # Define the pulses
        if 'exc_pulse' in self.user_inputs.keys():
            if self.user_inputs['exc_pulse'] == 'Rectangular':
                exc_pulse = RectPulse(  
                    tp=16, freq=0, flipangle=np.pi/2, scale=0
                )
            elif self.user_inputs['exc_pulse'] == 'Chirp':
                exc_pulse = ChirpPulse(
                    tp=16, init_freq=-0.25, final_freq=-0.03, flipangle=np.pi/2, scale=0,
                )
            elif self.user_inputs['exc_pulse'] == 'HS':
                exc_pulse = HSPulse(  
                    tp=120, init_freq=-0.25, final_freq=-0.03, flipangle=np.pi/2, scale=0,
                    order1=6, order2=1, beta=10
                )
        else:
            exc_pulse = RectPulse(  
                    tp=16, freq=0, flipangle=np.pi/2, scale=0
                )
            
        if 'ref_pulse' in self.user_inputs.keys():
            if self.user_inputs['ref_pulse'] == 'Rectangular':
                ref_pulse = RectPulse(  
                    tp=16, freq=0, flipangle=np.pi, scale=0
                )
            elif self.user_inputs['ref_pulse'] == 'Chirp':
                ref_pulse = ChirpPulse(
                    tp=16, init_freq=-0.25, final_freq=-0.03, flipangle=np.pi, scale=0,
                )
            elif self.user_inputs['ref_pulse'] == 'HS':
                ref_pulse = HSPulse(  
                    tp=120, init_freq=-0.25, final_freq=-0.03, flipangle=np.pi, scale=0,
                    order1=6, order2=1, beta=10
                )
        else:
            ref_pulse = RectPulse(  
                            tp=16, freq=0, flipangle=np.pi, scale=0
                        )
            
        if 'pump_pulse' in self.user_inputs.keys():
            if self.user_inputs['pump_pulse'] == 'Rectangular':
                pump_pulse = RectPulse(  
                    tp=16, freq=0, flipangle=np.pi, scale=0
                )
            elif self.user_inputs['pump_pulse'] == 'Chirp':
                pump_pulse = ChirpPulse(
                    tp=16, init_freq=-0.25, final_freq=-0.03, flipangle=np.pi, scale=0,
                )
            elif self.user_inputs['pump_pulse'] == 'HS':
                pump_pulse = HSPulse(  
                    tp=120, init_freq=-0.25, final_freq=-0.03, flipangle=np.pi, scale=0,
                    order1=6, order2=1, beta=10
                )
        else:
            pump_pulse = HSPulse(  
                        tp=120, init_freq=-0.25, final_freq=-0.03, flipangle=np.pi, scale=0,
                        order1=6, order2=1, beta=10
                    )
            
        det_event = Detection(tp=512, freq=0)

        self.pulses = {'exc_pulse':exc_pulse, 'ref_pulse':ref_pulse, 'pump_pulse':pump_pulse, 'det_event':det_event}

        # Optimise the pulse positions
        self.signals.status.emit('Optimising pulses')
        self.signals.optimise_pulses.emit(self.pulses)
        self.pause_and_wait()
        # Tune the pulses
        self.signals.status.emit('Tuning pulses')
        exc_pulse = self.interface.tune_pulse(exc_pulse, mode="amp_nut", B=LO/gyro_N,LO=LO,reptime=reptime,shots=100)
        ref_pulse = self.interface.tune_pulse(ref_pulse, mode="amp_nut", B=LO/gyro_N,LO=LO,reptime=reptime,shots=100)
        pump_pulse = self.interface.tune_pulse(pump_pulse, mode="amp_nut", B=LO/gyro_N,LO=LO,reptime=reptime,shots=100)

        # Run a relaxation experiment
        self.signals.status.emit('Running relaxation Experiment')
        self.run_relax(reptime)
        self.signals.status.emit('Analysing relaxation Experiment')

        self.pause_and_wait()

        if ('tau1' in self.user_inputs) and (self.user_inputs['tau1'] != 0):
            self.signals.status.emit('Skipping QuickDEER')
        else:
            self.signals.status.emit('Running QuickDEER')

            tau2hrs = self.results['relax'].tau2hrs

            self.run_quick_deer(reptime,tau=tau2hrs)
            self.signals.status.emit('Analysing QuickDEER')
            self.pause_and_wait()

        self.signals.status.emit('Running LongDEER')
        if ('tau1' in self.user_inputs) and (self.user_inputs['tau1'] != 0):
            tau1 = self.user_inputs['tau1']
            tau2 = self.user_inputs['tau2']
            tau3 = self.user_inputs['tau3']
            deertype = self.user_inputs['ExpType']
            self.run_long_deer(reptime,tau1,tau2,tau3,deertype =deertype)
        else:
            rec_tau = self.results['quickdeer'].rec_tau_max
            max_tau = self.results['relax'].max_tau
            tau = np.min([rec_tau,max_tau])
            self.run_long_deer(reptime,tau,tau,0.3,deertype ='5pDEER')

        self.signals.finished.emit()

    def new_data(self, data):
        self.results = data


