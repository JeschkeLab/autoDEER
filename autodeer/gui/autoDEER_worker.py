import PyQt6.QtCore as QtCore
from autodeer import DEERCriteria
from pyepr.sequences import *
from autodeer.sequences import *
from pyepr.criteria import *
from pyepr import get_waveform_precision
import pyepr as epr
import time
import numpy as np
from threadpoolctl import threadpool_limits
from collections import deque
import copy


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

    #Signal into autoDEER worker
    update_deer_settings = QtCore.pyqtSignal(object)
    

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

    def __init__(self, interface, wait:QtCore.QWaitCondition, 
                 mutex:QtCore.QMutex,pulses:dict, results:dict, freq, gyro, 
                 AWG=True, user_inputs:dict = None, 
                 operating_mode = None,fixed_tau=None, *args, **kwargs):
        super(autoDEERWorker,self).__init__()

        self.methods = deque()
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
        self.freq = freq
        self.gyro = gyro
        self.AWG = AWG
        self.user_inputs = user_inputs
        self.stop_flag = False
        self.quick_deer_state = True
        self.fixed_tau = fixed_tau
        self.operating_mode = operating_mode

        self.max_tau = 3.5
        print(f"Waveform Precision is: {get_waveform_precision()}")

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

        self.Q = None
        self.skip_list = []
        # # Add the callback to our kwargs
        # self.kwargs['progress_callback'] = self.signals.progress

        self.signals.update_deer_settings.connect(self.update_deersettings)

    def pause_and_wait(self):
        self.mutex.lock()
        self.wait.wait(self.mutex)
        self.mutex.unlock()


    def run_fsweep(self):
        '''
        Initialise the runner function with passed args, kwargs.
        '''
        freq = self.freq
        gyro_N = self.gyro
        reptime = self.reptime
        p90, p180 = self.interface.tune_rectpulse(tp=self.tp, freq=freq, B=freq/gyro_N, reptime = reptime,shots=int(100*self.noise_mode))
        shots = int(150*self.noise_mode)
        shots = np.max([shots,50])
        fsweep = FieldSweepSequence(
            B=freq/gyro_N, freq=freq,reptime=reptime,averages=50,shots=shots,
            Bwidth = 250, 
            pi2_pulse=p90, pi_pulse=p180,
        )

        self.interface.launch(fsweep,savename=self.savename("EDFS_Q"),)
        self.signals.status.emit('Field-sweep running')
        self.interface.terminate_at(SNRCriteria(150),test_interval=self.test_interval,verbosity=2)
        while self.interface.isrunning():
            time.sleep(self.updaterate)
        self.signals.status.emit('Field-sweep complete')
        self.signals.fsweep_result.emit(self.interface.acquire_dataset())

    def run_respro(self):
        '''
        Initialise the runner function for resonator profile.
        '''
        freq = self.freq
        gyro=self.gyro
        reptime = self.reptime
        p90, p180 = self.interface.tune_rectpulse(tp=self.tp, freq=freq, B=freq/gyro, reptime = reptime,shots=int(100*self.noise_mode))
        shots = int(200*self.noise_mode)
        shots = np.max([shots,25])
        dtp = np.max([1, int(np.round(get_waveform_precision()))])

        if self.Q is None or self.Q == 0:
            fwidth=0.3
        else:
            fwidth=np.around(freq/self.Q,2)
            
        RPseq = ResonatorProfileSequence(
            B=freq/gyro, freq=freq,reptime=reptime,averages=10,shots=shots,
            pi2_pulse=p90, pi_pulse=p180, fwidth=fwidth, dtp=dtp,
        )

        self.interface.launch(RPseq,savename=self.savename("ResPro"),)
        self.signals.status.emit('Resonator Profile running')

        self.interface.terminate_at(SNRCriteria(5),test_interval=self.test_interval,verbosity=2)
        while self.interface.isrunning():
            time.sleep(self.updaterate)
        self.signals.status.emit('Resonator Profile complete')
        self.signals.respro_result.emit(self.interface.acquire_dataset())
        # self.pause_and_wait()
        # if np.abs(freq-self.freq) > 0.1:
        #     # Rerun Resonator Profile
        #     self.run_respro()

        # return 'skip'
            


    def run_CP_relax(self,dt=20,tmin=0.5,averages=30,autoStop=True,autoIFGain=True,*kwargs):
        '''
        Initialise the runner function for relaxation. 
        '''
        self.signals.status.emit('Running Carr-Purcell Experiment')
        freq = self.freq
        gyro = self.gyro
        reptime = self.reptime
        shots = int(50*self.noise_mode)
        shots = np.max([shots,5])
        relax = DEERSequence(
            B=freq/gyro, freq=freq,reptime=reptime,averages=averages,shots=shots,
            tau1=tmin, tau2=tmin, tau3=0.3, dt=dt,
            exc_pulse=self.pulses['exc_pulse'], ref_pulse=self.pulses['ref_pulse'],
            pump_pulse=self.pulses['pump_pulse'], det_event=self.pulses['det_event']
            )
        print(f"Running Carr-Purcell Experiment with tmin: {relax.tau2.value*2} us and dt: {relax.dt} ns")

        relax.five_pulse(relaxation=True, re_step=dt, re_dim=250)
        if self.AWG:
            relax.select_pcyc("16step_5p")
        else:
            relax.select_pcyc("DC")
            relax.shots.value *= 8 
        relax._estimate_time()

        # relax.pulses[1].scale.value = 0
        # relax.pulses[3].scale.value = 0
        if not autoIFGain:
            print(f"Using manual IF Gain. IF Gain: {self.interface.IFgain}")
            autoIFGain = self.interface.IFgain
        else:
            autoIFGain = None

        self.interface.launch(relax, savename=self.savename("CP"), IFgain=autoIFGain)
        if autoStop:
            self.interface.terminate_at(SNRCriteria(30, verbosity=2), test_interval=self.test_interval, verbosity=2)
        while self.interface.isrunning():
            time.sleep(self.updaterate)
        self.signals.relax_result.emit(self.interface.acquire_dataset())
        self.signals.status.emit('Carr-Purcell experiment complete')
        
    def run_T2_relax(self,dt=20,tmin=0.4,averages=30,autoStop=True,autoIFGain=True,*kwargs):
        self.signals.status.emit('Running T2 experiment')
        freq = self.freq
        gyro = self.gyro
        reptime = self.reptime
        shots = int(100*self.noise_mode)
        shots = np.max([shots,15])

        seq = T2RelaxationSequence(
            B=freq/gyro, freq=freq,reptime=reptime,averages=averages,shots=shots,
            start=tmin*1e3,step=dt,dim=250,pi2_pulse=self.pulses['exc_pulse'],
            pi_pulse=self.pulses['ref_pulse'], det_event=self.pulses['det_event'])
        
        if not autoIFGain:
            print(f"Using manual IF Gain. IF Gain: {self.interface.IFgain}")
            autoIFGain = self.interface.IFgain
        else:
            autoIFGain = None

        self.interface.launch(seq,savename=self.savename("T2_Q"),IFgain=autoIFGain)
        if autoStop:
            self.interface.terminate_at(SNRCriteria(30),test_interval=self.test_interval,verbosity=2)
        while self.interface.isrunning():
            time.sleep(self.updaterate)
        self.signals.T2_result.emit(self.interface.acquire_dataset())
        self.signals.status.emit('T2relax experiment complete')

    def run_2D_relax(self):
        self.signals.status.emit('Running 2D decoherence experiment')
        freq = self.freq
        gyro = self.gyro
        reptime = self.reptime

        tau = self.max_tau
        dim = 75
        seq = RefocusedEcho2DSequence(
            B=freq/gyro, freq=freq,reptime=reptime,averages=10,shots=int(50*self.noise_mode),
            tau=tau,dim=dim, pi2_pulse=self.pulses['exc_pulse'],
            pi_pulse=self.pulses['ref_pulse'], det_event=self.pulses['det_event'])
        

        self.interface.launch(seq,savename=self.savename("2D_DEC"),)
        self.interface.terminate_at(SNRCriteria(15),test_interval=self.test_interval,verbosity=2,)
        while self.interface.isrunning():
            time.sleep(self.updaterate)
        self.signals.Relax2D_result.emit(self.interface.acquire_dataset())
        self.signals.status.emit('2D decoherence experiment complete')

    def run_1D_refocused_echo(self,dt=20,tmin=0.4,averages=30,autoStop=True,autoIFGain=True,*kwargs):

        self.signals.status.emit('Running 1D refocused echo experiment')
        freq = self.freq
        gyro = self.gyro
        reptime = self.reptime
        shots = int(25*self.noise_mode)
        shots = np.max([shots,5])

        seq = RefocusedEcho1DSequence(
            B=freq/gyro, freq=freq,reptime=reptime,averages=averages,shots=shots,
            tau1=400, start=tmin*1e3, step=dt, dim=250,
            pi2_pulse=self.pulses['exc_pulse'], pi_pulse=self.pulses['ref_pulse'], det_event=self.pulses['det_event'], pump_pulse=self.pulses['pump_pulse']
        )

        if not autoIFGain:
            autoIFGain = self.interface.IFgain
        else:
            autoIFGain = None

        self.interface.launch(seq,savename=self.savename("1DRefEcho"),IFgain=autoIFGain)

        if autoStop:
            self.interface.terminate_at(SNRCriteria(30),test_interval=self.test_interval,verbosity=2)
        while self.interface.isrunning():
            time.sleep(self.updaterate)
        self.signals.T2_result.emit(self.interface.acquire_dataset())
        self.signals.status.emit('1D refocused echo experiment experiment complete')

    def run_quick_deer(self):

        if not self.quick_deer_state:
            return 'skip'
        
        self.signals.status.emit('Running QuickDEER')
        criteria = self.deer_inputs['criteria']
        DEER_crit = DEERCriteria(mode=criteria,verbosity=2,update_func=self.signals.quickdeer_update.emit)
        time_crit = TimeCriteria('Max time criteria', time.time() + 4*60*60 ,'Max time criteria')
        total_crit = [DEER_crit, time_crit]
        signal = self.signals.quickdeer_result.emit
        self.run_deer(total_crit, signal, dt=16,shot=15,averages=150 )

    def run_long_deer(self):
        self.signals.status.emit('Running LongDEER')
        if 'autoStop' in self.deer_inputs and not self.deer_inputs['autoStop']:
            DEER_crit = DEERCriteria(mode=np.inf,verbosity=2,update_func=self.signals.longdeer_update.emit)
            total_crit = [DEER_crit]
        else: # autoStop is True
            DEER_crit = DEERCriteria(mode="high",verbosity=2,update_func=self.signals.longdeer_update.emit)
            total_crit = [DEER_crit, self.EndTimeCriteria]
        end_signal = self.signals.longdeer_result.emit
        self.run_deer(total_crit,end_signal, dt=16,shot=50,averages=1000)

    def run_single_deer(self):
        # Run a DEER experiment background measurement
        self.signals.status.emit('Running DEER Background')
        if 'autoStop' in self.deer_inputs and not self.deer_inputs['autoStop']:
            SNR_crit = SNRCriteria(150,verbosity=2)
            total_crit = [SNR_crit]
        else: # autoStop is True
            SNR_crit = SNRCriteria(150,verbosity=2)
            total_crit = [SNR_crit, self.EndTimeCriteria]
        end_signal = self.signals.longdeer_result.emit
        self.run_deer(total_crit,end_signal, dt=16,shot=50,averages=1000)

    def run_deer(self, end_criteria, signal, dt=16, shot=50, averages=1000,):
    
        freq = self.freq
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
            tau = np.min([rec_tau, max_tau])
            deertype = '5pDEER'
            tau1 = tau
            tau2 = tau
            tau3 = 0.3

        if 'ESEEM' in self.deer_inputs:
            ESEEM = self.deer_inputs['ESEEM']
        else:
            ESEEM = None

        deer = DEERSequence(
            B=freq/self.gyro, freq=freq, reptime=reptime, averages=averages,
            shots=int(250*self.noise_mode), tau1=tau1, tau2=tau2, tau3=tau3,
            dt=dt, exc_pulse=self.pulses['exc_pulse'],
            ref_pulse=self.pulses['ref_pulse'],
            pump_pulse=self.pulses['pump_pulse'],
            det_event=self.pulses['det_event'],
            ESEEM_avg=ESEEM
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

        self.interface.launch(deer, savename=self.savename(savename_type, savename_suffix),)

        time.sleep(30)  # Always wait for the experiment to properly start
        with threadpool_limits(limits=self.cores, user_api='blas'):
            self.interface.terminate_at(end_criteria, verbosity=2, test_interval=self.test_interval) # Change criteria backagain
        while self.interface.isrunning():
            time.sleep(self.updaterate)
        signal(self.interface.acquire_dataset())
        self.signals.status.emit('DEER experiment complete')

    def run_reptime_opt(self):
        reptime_guess = self.reptime
        freq = self.freq
        p90, p180 = self.interface.tune_rectpulse(tp=self.tp, freq=freq, B=freq/self.gyro, reptime = reptime_guess,shots=int(100*self.noise_mode))

        n_shots = int(np.max([int(50*self.noise_mode),10]))
        scan = ReptimeScan(B=freq/self.gyro, freq=freq,reptime=reptime_guess, reptime_max=12e3, averages=10, shots=n_shots,
                           pi2_pulse=p90, pi_pulse=p180)
        self.interface.launch(scan,savename=f"{self.samplename}_reptimescan",)
        self.interface.terminate_at(SNRCriteria(30),verbosity=2,test_interval=self.test_interval)
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
        shots = np.max([int(10*self.noise_mode), 2])
        self.signals.status.emit('Tuning pulses')
        exc_pulse = self.interface.tune_pulse(exc_pulse, mode="amp_nut", B=self.freq/self.gyro,freq=self.freq,reptime=self.reptime,shots=shots)
        ref_pulse = self.interface.tune_pulse(ref_pulse, mode="amp_nut", B=self.freq/self.gyro,freq=self.freq,reptime=self.reptime,shots=shots)
        if isinstance(pump_pulse, epr.FrequencySweptPulse):  # A frequency swept pump pulse's optinmal power is max
            amp_factor = pump_pulse.amp_factor.value
            freq_range = np.linspace(pump_pulse.init_freq.value, pump_pulse.final_freq.value, 100) + self.freq
            B1 = self.interface.resonator.model_func(freq_range)

            scale = amp_factor/np.max(B1)
            if scale > 1:
                scale = 1
            pump_pulse.scale.value = scale
        else:
            pump_pulse = self.interface.tune_pulse(pump_pulse, mode="amp_nut", B=self.freq/self.gyro,freq=self.freq,reptime=self.reptime,shots=shots)

        return 'skip'

    def _build_methods(self):
        
        methods = deque([self.run_fsweep,self.run_reptime_opt,self.run_respro])
        
        if ( self.operating_mode is None) or (self.operating_mode == '5pDEER') or (self.operating_mode == 'auto'):
            methods.append(self.run_CP_relax)
            methods.append(self.run_1D_refocused_echo)
            methods.append(self.run_T2_relax)
        
        elif ( self.operating_mode == '4pDEER') or ( self.operating_mode == 'Ref2D'):
            methods.append(self.run_CP_relax)
            methods.append(self.run_1D_refocused_echo)
            methods.append(self.run_T2_relax)
            methods.append(self.run_2D_relax)
            if  self.operating_mode == 'Ref2D':
                return methods

        elif  self.operating_mode == 'single':
            methods.append(self.run_T2_relax)
            methods.append(self.run_single_deer)
            return methods        

        if self.fixed_tau is None:
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
        
        self.methods = self._build_methods()
        if self.pulses != {} or ('exc_pulses' in self.pulses and self.pulses['exc_pulse'].scale.value == 0):
            self.methods.appendleft(self.tune_pulses)

        print(f'Skip List:', self.skip_list)
        while self.methods:
            if self.stop_flag:
                self.signals.finished.emit()
                return None

            method = self.methods.popleft()
            if method.__name__ in self.skip_list:
                print(f"Skipping",method.__name__)
                self.skip_list.remove(method.__name__)
                continue

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
        print("New pulses in autoDEER worker")
        self.pulses = pulses
        self.methods.appendleft(self.tune_pulses)


    def update_reptime(self,reptime):
        self.reptime = reptime

    def update_gyro(self,gyro):
        self.gyro = gyro

    def set_2D_max_tau(self, max_tau):
        self.max_tau = max_tau

    def update_deersettings(self,deer_settings):
        self.deer_inputs = copy.deepcopy(deer_settings)

    def stop(self):
        self.stop_flag = True
        self.methods.clear()
        self.signals.status.emit('Stopping experiment')
        
    def set_noise_mode(self,level):
        self.noise_mode = level
    
    def update_freq(self,freq,repeat=False):
        self.freq = freq

        if repeat:
            self.methods.appendleft(self.run_respro)


    def repeat_fieldsweep(self):
        self.methods.appendleft(self.run_fsweep)

    def repeat_quickdeer(self):
        self.methods.appendleft(self.run_quick_deer)

    def rerun_relax(self,experiment,**kwargs):
        if experiment == 'CP-relax':
            method = lambda: self.run_CP_relax(**kwargs)
        elif experiment == 'Tm-relax':
            method = lambda: self.run_T2_relax(**kwargs)
        elif experiment == 'RefEcho2D':
            method = lambda: self.run_2D_relax(**kwargs)
        elif experiment == 'RefEcho1D':
            method = lambda: self.run_1D_refocused_echo(**kwargs)
        self.methods.appendleft(method)

    def setQ(self,Q):
        self.Q = Q

    def update_skip_list(self,skip_list):
        self.skip_list = skip_list
        # old_methods = self.methods

        # new_methods = deque(item for item in old_methods if not (hasattr(item,'__name__') and item.__name__ in skip_list))
        # self.methods = new_methods
        # print('Old methods', list(old_methods))
        # print('New methods', list(new_methods))