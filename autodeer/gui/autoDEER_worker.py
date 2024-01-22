import PyQt6.QtCore as QtCore
from autodeer import RectPulse, ChirpPulse, HSPulse, Detection, DEERCriteria, SNRCriteria, TimeCriteria
from autodeer.sequences import *
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
    # optimise_pulses = QtCore.pyqtSignal(object)
    relax_result = QtCore.pyqtSignal(object)
    T2_result = QtCore.pyqtSignal(object)
    quickdeer_result = QtCore.pyqtSignal(object)
    longdeer_result = QtCore.pyqtSignal(object)
    quickdeer_update = QtCore.pyqtSignal(object)
    longdeer_update = QtCore.pyqtSignal(object)
    reptime_scan_result = QtCore.pyqtSignal(object)
    timeout = QtCore.pyqtSignal()
    

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

    def __init__(self, interface, wait:QtCore.QWaitCondition, mutex:QtCore.QMutex,pulses:dict, results:dict, LO, gyro,AWG=True, user_inputs:dict = None, *args, **kwargs):
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
        self.pulses = pulses
        self.samplename = user_inputs['sample']
        self.LO = LO
        self.gyro = gyro
        self.AWG = AWG
        self.user_inputs = user_inputs
        self.stop_flag = False

        if 'SampleConc' in self.user_inputs:
            self.noise_mode = self.user_inputs['SampleConc']
        else:
            self.noise_mode = 1
        

        if not 'DEER_update_func' in self.user_inputs:
            self.user_inputs['DEER_update_func'] = None
        
        if 'cores' in kwargs:
            self.cores = kwargs['cores']
        else:
            self.cores = 1

        self.deer_inputs = {}
        
        self.EndTimeCriteria = TimeCriteria('End Time',time.time() + self.user_inputs['MaxTime']*3600, "Overall end time",end_signal=self.signals.timeout.emit)

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
        LO = self.LO
        gyro_N = self.gyro
        reptime = self.reptime
        p90, p180 = self.interface.tune_rectpulse(tp=12, LO=LO, B=LO/gyro_N, reptime = reptime,shots=int(100*self.noise_mode))
        fsweep = FieldSweepSequence(
            B=LO/gyro_N, LO=LO,reptime=reptime,averages=10,shots=int(150*self.noise_mode),
            Bwidth = 350, 
            pi2_pulse=p90, pi_pulse=p180,
        )

        self.interface.launch(fsweep,savename=f"{self.samplename}_fieldsweep",IFgain=2)
        self.signals.status.emit('Field-sweep running')
        self.interface.terminate_at(SNRCriteria(10))
        # while self.interface.isrunning():
        #     time.sleep(self.updaterate)
        self.signals.status.emit('Field-sweep complete')
        self.signals.fsweep_result.emit(self.interface.acquire_dataset())

    def run_respro(self):
        '''
        Initialise the runner function for resonator profile.
        '''
        LO = self.LO
        gyro=self.gyro
        reptime = self.reptime
        p90, p180 = self.interface.tune_rectpulse(tp=12, LO=LO, B=LO/gyro, reptime = reptime,shots=int(100*self.noise_mode))

        RPseq = ResonatorProfileSequence(
            B=LO/gyro, LO=LO,reptime=reptime,averages=10,shots=int(40*self.noise_mode),
            pi2_pulse=p90, pi_pulse=p180,
        )

        self.interface.launch(RPseq,savename=f"{self.samplename}_resonator_profile",IFgain=2)
        self.signals.status.emit('Resonator Profile running')

        self.interface.terminate_at(SNRCriteria(5))
        # while self.interface.isrunning():
        #     time.sleep(self.updaterate)
        self.signals.status.emit('Resonator Profile complete')
        self.signals.respro_result.emit(self.interface.acquire_dataset())
        self.pause_and_wait()
        if np.abs(LO-self.LO) > 0.1:
            # Rerun Resonator Profile
            self.run_respro()

        return 'skip'
            


    def run_CP_relax(self):
        '''
        Initialise the runner function for relaxation. 
        '''
        self.signals.status.emit('Running Carr-Purcell Experiment')
        LO = self.LO
        gyro = self.gyro
        reptime = self.reptime
        relax = DEERSequence(
            B=LO/gyro, LO=LO,reptime=reptime,averages=4,shots=int(50*self.noise_mode),
            tau1=0.5, tau2=0.5, tau3=0.2, dt=15,
            exc_pulse=self.pulses['exc_pulse'], ref_pulse=self.pulses['ref_pulse'],
            pump_pulse=self.pulses['pump_pulse'], det_event=self.pulses['det_event']
            )
        relax.five_pulse(relaxation=True, re_step=200)
        if self.AWG:
            relax.select_pcyc("16step_5p")
        else:
            relax.select_pcyc("DC")
        relax._estimate_time();
        relax.pulses[1].scale.value = 0
        relax.pulses[3].scale.value = 0
        self.interface.launch(relax,savename=f"{self.samplename}_CP2relax",IFgain=2)
        # self.interface.terminate_at(SNRCriteria(30),test_interval=0.5)
        while self.interface.isrunning():
            time.sleep(self.updaterate)
        self.signals.relax_result.emit(self.interface.acquire_dataset())
        self.signals.status.emit('Carr-Purcell experiment complete')
        
    def run_T2_relax(self):
        self.signals.status.emit('Running T2 experiment')
        LO = self.LO
        gyro = self.gyro
        reptime = self.reptime

        seq = T2RelaxationSequence(
            B=LO/gyro, LO=LO,reptime=reptime,averages=4,shots=int(50*self.noise_mode),
            step=60,dim=200,pi2_pulse=self.pulses['exc_pulse'],
            pi_pulse=self.pulses['ref_pulse'], det_event=self.pulses['det_event'])
        
        self.interface.launch(seq,savename=f"{self.samplename}_T2relax",IFgain=2)
        while self.interface.isrunning():
            time.sleep(self.updaterate)
        self.signals.T2_result.emit(self.interface.acquire_dataset())
        self.signals.status.emit('T2relax experiment complete')


    def run_2D_relax(self):
        self.signals.status.emit('Running 2D relaxation experiment')
        LO = self.LO
        gyro = self.gyro
        reptime = self.reptime

    def run_quick_deer(self, dt=16):
        if ('tau1' in self.user_inputs) and (self.user_inputs['tau1'] != 0):
            self.signals.status.emit('Skipping QuickDEER')
            return 'skip'
        else:
            self.signals.status.emit('Running QuickDEER')
            
            tau = self.results['relax'].tau2hrs
            
            LO = self.LO
            gyro = self.gyro
            reptime = self.reptime
            deer = DEERSequence(
                B=LO/gyro, LO=LO,reptime=reptime,averages=150,shots=int(15*self.noise_mode),
                tau1=tau, tau2=tau, tau3=0.3, dt=dt,
                exc_pulse=self.pulses['exc_pulse'], ref_pulse=self.pulses['ref_pulse'],
                pump_pulse=self.pulses['pump_pulse'], det_event=self.pulses['det_event']
                )
            
            deer.five_pulse()
            if self.AWG:
                deer.select_pcyc("16step_5p")
            else:
                deer.select_pcyc("DC")
            deer._estimate_time();

            # build criteria
            DEER_crit = DEERCriteria(mode="speed",verbosity=2,update_func=self.signals.quickdeer_update.emit)
            time_crit = TimeCriteria('Max time criteria', time.time() + 4*60*60 ,'Max time criteria')
            total_crit = DEER_crit + time_crit

            self.interface.launch(deer,savename=f"{self.samplename}_quickdeer",IFgain=2)
            with threadpool_limits(limits=self.cores, user_api='blas'):
                self.interface.terminate_at(total_crit,verbosity=2,test_interval=0.5)
            self.signals.quickdeer_result.emit(self.interface.acquire_dataset())
            self.signals.status.emit('QuickDEER experiment complete')

    def run_long_deer(self, dt=16, deertype='5pDEER'):
        self.signals.status.emit('Running LongDEER')
        LO = self.LO
        reptime = self.reptime

        if ('tau1' in self.user_inputs) and (self.user_inputs['tau1'] != 0):
            tau1 = self.user_inputs['tau1']
            tau2 = self.user_inputs['tau2']
            tau3 = self.user_inputs['tau3']
            deertype = self.user_inputs['ExpType']
        elif self.deer_inputs != {}:
            tau1 = self.deer_inputs['tau1']
            tau2 = self.deer_inputs['tau2']
            tau3 = self.deer_inputs['tau3']
            deertype = self.deer_inputs['ExpType']
            dt = self.deer_inputs['dt']
 
        else:
            rec_tau = self.results['quickdeer'].rec_tau_max
            dt = self.results['quickdeer'].rec_dt * 1e3
            max_tau = self.results['relax'].max_tau
            tau = np.min([rec_tau,max_tau])
            tau1 = tau
            tau2 = tau
            tau3 = 0.3


        deer = DEERSequence(
            B=LO/self.gyro, LO=LO,reptime=reptime,averages=1000,shots=int(50*self.noise_mode),
            tau1=tau1, tau2=tau2, tau3=tau3, dt=dt,
            exc_pulse=self.pulses['exc_pulse'], ref_pulse=self.pulses['ref_pulse'],
            pump_pulse=self.pulses['pump_pulse'], det_event=self.pulses['det_event'],
            ESEEM_avg = self.deer_inputs['ESEEM']
            )
        
        if deertype == '5pDEER':
            deer.five_pulse()
            savename_suffix = f"5p_{tau1:.3f}us_{tau2:.3f}us_{tau3:.3f}us"
        elif deertype == '4pDEER':
            deer.four_pulse()
            savename_suffix = f"4p_{tau1:.3f}us_{tau2:.3f}us"
        elif deertype == '3pDEER':
            deer.four_pulse()
            savename_suffix = f"3p_{tau1:.3f}u"
        
        if not self.AWG:
            deer.select_pcyc('DC')
        elif deertype == '5pDEER':
            deer.select_pcyc("16step_5p")
        elif deertype == '4pDEER':
            deer.select_pcyc("16step_4p")
        elif deertype == '3pDEER':
            deer.select_pcyc("8step_3p")

        deer._estimate_time();
        # build criteria
        DEER_crit = DEERCriteria(mode="high",verbosity=2,update_func=self.signals.longdeer_update.emit)
        total_crit = DEER_crit + self.EndTimeCriteria

        self.interface.launch(deer,savename=f"{self.samplename}_deer"+savename_suffix,IFgain=2)
        time.sleep(30) # Always wait for the experiment to properly start
        with threadpool_limits(limits=self.cores, user_api='blas'):
            self.interface.terminate_at(total_crit,verbosity=2,test_interval=0.5) # Change criteria backagain
        self.signals.longdeer_result.emit(self.interface.acquire_dataset())
        self.signals.status.emit('Long DEER experiment complete')

    def run_reptime_opt(self):
        reptime_guess = self.reptime
        LO = self.LO
        p90, p180 = self.interface.tune_rectpulse(tp=12, LO=LO, B=LO/self.gyro, reptime = reptime_guess,shots=int(100*self.noise_mode))

        scan = ReptimeScan(B=LO/self.gyro, LO=LO,reptime=reptime_guess, reptime_max=12e3, averages=10, shots=int(50*self.noise_mode),
                           pi2_pulse=p90, pi_pulse=p180)
        self.interface.launch(scan,savename=f"{self.samplename}_reptimescan",IFgain=2)
        self.interface.terminate_at(SNRCriteria(15),verbosity=2,test_interval=0.5)
        # while self.interface.isrunning():
        #     time.sleep(self.updaterate)
        self.signals.status.emit('Reptime scan complete')
        self.signals.reptime_scan_result.emit(self.interface.acquire_dataset())

    def tune_pulses(self):
        # Tune the pulses
        self.signals.status.emit('Tuning pulses')
        pump_pulse = self.pulses['pump_pulse']
        ref_pulse = self.pulses['ref_pulse']
        exc_pulse = self.pulses['exc_pulse']
        det_event = self.pulses['det_event']
        
        self.signals.status.emit('Tuning pulses')
        exc_pulse = self.interface.tune_pulse(exc_pulse, mode="amp_nut", B=self.LO/self.gyro,LO=self.LO,reptime=self.reptime,shots=int(100*self.noise_mode))
        ref_pulse = self.interface.tune_pulse(ref_pulse, mode="amp_nut", B=self.LO/self.gyro,LO=self.LO,reptime=self.reptime,shots=int(100*self.noise_mode))
        pump_pulse = self.interface.tune_pulse(pump_pulse, mode="amp_nut", B=self.LO/self.gyro,LO=self.LO,reptime=self.reptime,shots=int(100*self.noise_mode))

        return 'skip'

    @QtCore.pyqtSlot()    
    def run(self):

        self.reptime = 3e3
        self.stop_flag = False


        methods = [self.run_fsweep,self.run_reptime_opt,self.run_respro,self.run_fsweep,
                   self.tune_pulses,self.run_CP_relax,self.run_T2_relax,self.run_quick_deer,
                   self.run_long_deer]

        for method in methods:
            if self.stop_flag:
                self.signals.finished.emit()
                return None
            flag = method()
            if flag is None:
                self.pause_and_wait()
            elif flag == 'skip':
                continue
            
            if self.stop_flag:
                self.signals.finished.emit()
                return None
            
        self.signals.finished.emit()


    def new_data(self, data):
        self.results = data

    def new_pulses(self, pulses):
        self.pulses = pulses

    def update_LO(self,LO):
        self.LO = LO

    def update_reptime(self,reptime):
        self.reptime = reptime

    def update_deersettings(self,deer_settings):
        self.deer_inputs = deer_settings

    def stop(self):
        self.stop_flag = True
        self.signals.status.emit('Stopping experiment')
        

