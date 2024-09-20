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
    Relax2D_result = QtCore.pyqtSignal(object)
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
        self.project = user_inputs['project']
        self.LO = LO
        self.gyro = gyro
        self.AWG = AWG
        self.user_inputs = user_inputs
        self.stop_flag = False
        self.quick_deer_state = True

        self.max_tau = 3.5

        if (self.project is None) or (self.project == ''):
            def savename(exp, suffix=""):
                if suffix != "":
                    return f"({self.samplename})_({exp})_{suffix}"
                else:
                    return f"({self.samplename})_({exp})"
        else:
            def savename(exp,suffix=""):
                if suffix != "":
                    return f"({self.project})_({self.samplename})_({exp})_{suffix}"
                else:
                    return f"({self.project})_({self.samplename})_({exp})"
        self.savename = savename
        

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

        if 'tp' in kwargs:
            self.tp = kwargs['tp']
        else:
            self.tp=12

        if "night_hours" in kwargs:
            night_hours = kwargs['night_hours']

        self.deer_inputs = {}
        
        self.EndTimeCriteria = TimeCriteria('End Time',time.time() + self.user_inputs['MaxTime']*3600, "Overall end time",end_signal=self.signals.timeout.emit,night_hours=night_hours)
        self.test_interval = 0.5
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
        p90, p180 = self.interface.tune_rectpulse(tp=self.tp, LO=LO, B=LO/gyro_N, reptime = reptime,shots=int(100*self.noise_mode))
        shots = int(50*self.noise_mode)
        shots = np.min([shots,20])
        fsweep = FieldSweepSequence(
            B=LO/gyro_N, LO=LO,reptime=reptime,averages=50,shots=int(50*self.noise_mode),
            Bwidth = 250, 
            pi2_pulse=p90, pi_pulse=p180,
        )

        self.interface.launch(fsweep,savename=self.savename("EDFS_Q"),IFgain=1)
        self.signals.status.emit('Field-sweep running')
        self.interface.terminate_at(SNRCriteria(30),test_interval=self.test_interval)
        while self.interface.isrunning():
            time.sleep(self.updaterate)
        self.signals.status.emit('Field-sweep complete')
        self.signals.fsweep_result.emit(self.interface.acquire_dataset())

    def run_respro(self):
        '''
        Initialise the runner function for resonator profile.
        '''
        LO = self.LO
        gyro=self.gyro
        reptime = self.reptime
        p90, p180 = self.interface.tune_rectpulse(tp=self.tp, LO=LO, B=LO/gyro, reptime = reptime,shots=int(100*self.noise_mode))

        RPseq = ResonatorProfileSequence(
            B=LO/gyro, LO=LO,reptime=reptime,averages=10,shots=int(50*self.noise_mode),
            pi2_pulse=p90, pi_pulse=p180,fwidth=0.15
        )

        self.interface.launch(RPseq,savename=self.savename("ResPro"),IFgain=1)
        self.signals.status.emit('Resonator Profile running')

        self.interface.terminate_at(SNRCriteria(5),test_interval=self.test_interval)
        while self.interface.isrunning():
            time.sleep(self.updaterate)
        self.signals.status.emit('Resonator Profile complete')
        self.signals.respro_result.emit(self.interface.acquire_dataset())
        self.pause_and_wait()
        if np.abs(LO-self.LO) > 0.1:
            # Rerun Resonator Profile
            self.run_respro()

        return 'skip'
            


    def run_CP_relax(self,dt=200):
        '''
        Initialise the runner function for relaxation. 
        '''
        self.signals.status.emit('Running Carr-Purcell Experiment')
        LO = self.LO
        gyro = self.gyro
        reptime = self.reptime
        shots = int(40*self.noise_mode)
        shots = np.min([shots,10])
        relax = DEERSequence(
            B=LO/gyro, LO=LO,reptime=reptime,averages=10,shots=shots,
            tau1=0.5, tau2=0.5, tau3=0.2, dt=15,
            exc_pulse=self.pulses['exc_pulse'], ref_pulse=self.pulses['ref_pulse'],
            pump_pulse=self.pulses['pump_pulse'], det_event=self.pulses['det_event']
            )
        relax.five_pulse(relaxation=True, re_step=dt)
        if self.AWG:
            relax.select_pcyc("16step_5p")
        else:
            relax.select_pcyc("DC")
            relax.shots.value *= 8 
        relax._estimate_time();
        relax.pulses[1].scale.value = 0
        relax.pulses[3].scale.value = 0
        self.interface.terminate_at(SNRCriteria(50),test_interval=self.test_interval)
        self.interface.launch(relax,savename=self.savename("CP"),IFgain=1)
        while self.interface.isrunning():
            time.sleep(self.updaterate)
        self.signals.relax_result.emit(self.interface.acquire_dataset())
        self.signals.status.emit('Carr-Purcell experiment complete')
        
    def run_T2_relax(self,dt=60):
        self.signals.status.emit('Running T2 experiment')
        LO = self.LO
        gyro = self.gyro
        reptime = self.reptime
        shots = int(40*self.noise_mode)
        shots = np.min([shots,10])

        seq = T2RelaxationSequence(
            B=LO/gyro, LO=LO,reptime=reptime,averages=10,shots=shots,
            step=dt,dim=200,pi2_pulse=self.pulses['exc_pulse'],
            pi_pulse=self.pulses['ref_pulse'], det_event=self.pulses['det_event'])
        
        self.interface.launch(seq,savename=self.savename("T2_Q"),IFgain=1)
        self.interface.terminate_at(SNRCriteria(50),test_interval=0.5)
        while self.interface.isrunning():
            time.sleep(self.updaterate)
        self.signals.T2_result.emit(self.interface.acquire_dataset())
        self.signals.status.emit('T2relax experiment complete')


    def run_2D_relax(self):
        self.signals.status.emit('Running 2D decoherence experiment')
        LO = self.LO
        gyro = self.gyro
        reptime = self.reptime

        seq = RefocusedEcho2DSequence(
            B=LO/gyro, LO=LO,reptime=reptime,averages=10,shots=int(25*self.noise_mode),
            tau=self.max_tau, pi2_pulse=self.pulses['exc_pulse'],
            pi_pulse=self.pulses['ref_pulse'], det_event=self.pulses['det_event'])

        self.interface.launch(seq,savename=self.savename("2D_DEC"),IFgain=2)
        self.interface.terminate_at(SNRCriteria(15),test_interval=self.test_interval)
        while self.interface.isrunning():
            time.sleep(self.updaterate)
        self.signals.Relax2D_result.emit(self.interface.acquire_dataset())
        self.signals.status.emit('2D decoherence experiment complete')
    
    def run_quick_deer(self):

        if not self.quick_deer_state:
            return 'skip'
        
        self.signals.status.emit('Running QuickDEER')
        DEER_crit = DEERCriteria(mode="speed",verbosity=2,update_func=self.signals.quickdeer_update.emit)
        time_crit = TimeCriteria('Max time criteria', time.time() + 4*60*60 ,'Max time criteria')
        total_crit = [DEER_crit, time_crit]
        signal = self.signals.quickdeer_result.emit
        self.run_deer(total_crit,signal, dt=16,shot=15,averages=150 )

    def run_long_deer(self):
        self.signals.status.emit('Running LongDEER')
        DEER_crit = DEERCriteria(mode="high",verbosity=2,update_func=self.signals.longdeer_update.emit)
        total_crit = [DEER_crit, self.EndTimeCriteria]
        signal = self.signals.longdeer_result.emit
        self.run_deer(total_crit,signal, dt=16,shot=50,averages=1e4)


    def run_deer(self,end_criteria,signal, dt=16,shot=50,averages=1000,):
    
        LO = self.LO
        reptime = self.reptime

        if ('tau1' in self.user_inputs) and (self.user_inputs['tau1'] != 0):
            tau1 = self.user_inputs['tau1']
            tau2 = self.user_inputs['tau2']
            deertype = self.user_inputs['ExpType']
            if deertype == '5pDEER':
                tau3 = self.user_inputs['tau3']
            else:
                tau3 = None
            
        elif self.deer_inputs != {}:
            tau1 = self.deer_inputs['tau1']
            tau2 = self.deer_inputs['tau2']
            deertype = self.deer_inputs['ExpType']
            if deertype == '5pDEER':
                tau3 = self.deer_inputs['tau3']
            else:
                tau3 = None
            dt = self.deer_inputs['dt']
 
        else:
            rec_tau = self.results['quickdeer'].rec_tau_max
            dt = self.results['quickdeer'].rec_dt * 1e3
            max_tau = self.results['relax'].max_tau
            tau = np.min([rec_tau,max_tau])
            deertype = '5pDEER'
            tau1 = tau
            tau2 = tau
            tau3 = 0.3

        if 'ESEEM' in self.deer_inputs:
            ESEEM = self.deer_inputs['ESEEM']
        else:
            ESEEM = None


        deer = DEERSequence(
            B=LO/self.gyro, LO=LO,reptime=reptime,averages=averages,shots=int(50*self.noise_mode),
            tau1=tau1, tau2=tau2, tau3=tau3, dt=dt,
            exc_pulse=self.pulses['exc_pulse'], ref_pulse=self.pulses['ref_pulse'],
            pump_pulse=self.pulses['pump_pulse'], det_event=self.pulses['det_event'],
            ESEEM_avg = ESEEM
            )
        
        if deertype == '5pDEER':
            deer.five_pulse()
            savename_type = 'DEER_5P_Q'
            savename_suffix = f"{tau1:.3f}us_{tau2:.3f}us_{tau3:.3f}us"
        elif deertype == '4pDEER':
            deer.four_pulse()
            savename_type = 'DEER_4P_Q'
            savename_suffix = f"{tau1:.3f}us_{tau2:.3f}us"
        elif deertype == '3pDEER':
            deer.four_pulse()
            savename_type = 'DEER_3P_Q'
            savename_suffix = f"{tau1:.3f}u"
        
        if not self.AWG:
            deer.select_pcyc('DC')
            deer.shots.value *= 8 
        elif deertype == '5pDEER':
            deer.select_pcyc("16step_5p")
        elif deertype == '4pDEER':
            deer.select_pcyc("16step_4p")
        elif deertype == '3pDEER':
            deer.select_pcyc("8step_3p")

        deer._estimate_time();

        self.interface.launch(deer,savename=self.savename(savename_type,savename_suffix),IFgain=2)
        time.sleep(30) # Always wait for the experiment to properly start
        with threadpool_limits(limits=self.cores, user_api='blas'):
            self.interface.terminate_at(end_criteria,verbosity=2,test_interval=self.test_interval) # Change criteria backagain
        while self.interface.isrunning():
            time.sleep(self.updaterate)
        signal(self.interface.acquire_dataset())
        self.signals.status.emit('DEER experiment complete')

    def run_reptime_opt(self):
        reptime_guess = self.reptime
        LO = self.LO
        p90, p180 = self.interface.tune_rectpulse(tp=self.tp, LO=LO, B=LO/self.gyro, reptime = reptime_guess,shots=int(100*self.noise_mode))

        n_shots = int(np.min([int(50*self.noise_mode),50]))
        scan = ReptimeScan(B=LO/self.gyro, LO=LO,reptime=reptime_guess, reptime_max=12e3, averages=10, shots=n_shots,
                           pi2_pulse=p90, pi_pulse=p180)
        self.interface.launch(scan,savename=f"{self.samplename}_reptimescan",IFgain=1)
        self.interface.terminate_at(SNRCriteria(15),verbosity=2,test_interval=self.test_interval)
        while self.interface.isrunning():
            time.sleep(self.updaterate)
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

    def _build_methods(self):

        seq = self.deer_inputs['ExpType']

        if (self.deer_inputs['tau1'] == 0) and (self.deer_inputs['tau2'] == 0):
            quick_deer = True
        else:
            quick_deer = False
        
        methods = [self.run_fsweep,self.run_reptime_opt,self.run_respro,self.run_fsweep,
                   self.tune_pulses]
        
        if (seq is None) or (seq == '5pDEER'):
            methods.append(self.run_CP_relax)
            methods.append(self.run_T2_relax)
        elif (seq == '4pDEER') or (seq == 'Ref2D'):
            methods.append(self.run_CP_relax)
            methods.append(self.run_T2_relax)
            methods.append(self.run_2D_relax)

        if seq == 'Ref2D':
            return methods

        if quick_deer:
            methods.append(self.run_quick_deer)
        
        methods.append(self.run_long_deer)

        return methods

    @QtCore.pyqtSlot()    
    def run(self):

        self.reptime = 3e3
        self.stop_flag = False


        # methods = [self.run_fsweep,self.run_reptime_opt,self.run_respro,self.run_fsweep,
        #            self.tune_pulses,self.run_CP_relax,self.run_T2_relax,self.run_quick_deer,
        #            self.run_long_deer]
        
        methods = self._build_methods()

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

    def update_gyro(self,gyro):
        self.gyro = gyro

    def set_2D_max_tau(self, max_tau):
        self.max_tau = max_tau

    def update_deersettings(self,deer_settings):
        self.deer_inputs = deer_settings

    def stop(self):
        self.stop_flag = True
        self.signals.status.emit('Stopping experiment')
        
    def set_noise_mode(self,level):
        self.noise_mode = level

