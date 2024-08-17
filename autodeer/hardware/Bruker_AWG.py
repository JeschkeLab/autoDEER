from autodeer.classes import Interface, Parameter
from autodeer.pulses import Delay, Detection, RectPulse
from autodeer.sequences import Sequence, HahnEchoSequence
from autodeer.utils import save_file, transpose_list_of_dicts, transpose_dict_of_list, round_step
from autodeer import create_dataset_from_sequence, create_dataset_from_axes

from autodeer.hardware.XeprAPI_link import XeprAPILink
from autodeer.hardware.Bruker_tools import PulseSpel, run_general,build_unique_progtable,PSPhaseCycle, write_pulsespel_file
from autodeer.hardware.Bruker_MPFU import step_parameters

import tempfile
import time
from scipy.optimize import minimize_scalar, curve_fit
import numpy as np
import threading
import concurrent.futures
from PyQt6.QtCore import QThreadPool
from autodeer.gui import Worker
import os
from pathlib import Path
import datetime
from deerlab import correctphase
import matplotlib.pyplot as plt
import logging
import warnings
# ===================


# =============================================================================
hw_log = logging.getLogger('interface.Xepr')


class BrukerAWG(Interface):
    """
    Represents the interface for connecting to AWG based Bruker ELEXSYS-II 
    Spectrometers.
    """

    def __init__(self, config_file:str) -> None:
        """An interface for connecting to AWG based Bruker ELEXSYS-II 
        Spectrometers.

        Getting Started
        ------------------
        Before a connection can be made an appropriate configuration file first
        needs to be written. 

        1. Open Xepr
        2. Processing -> XeprAPI -> Enable XeprAPI
        3. `BrukerAWG.connect()`

        Parameters
        ----------
        config_file : str
            The path to a YAML configuration file.
        
        Attributes
        ----------
        bg_thread: None or threading.Thread
            If a background thread is needed, it is stored here.

        """

        self.api = XeprAPILink(config_file)
        self.spec_config = self.api.config["Spectrometer"]
        self.bridge_config = self.api.spec_config["Bridge"]
        
        self.temp_dir = tempfile.mkdtemp("autoDEER")

        self.d0 = self.bridge_config["d0"]
        
        self.bg_thread = None
        self.bg_data = None
        self.cur_exp = None
        self.tuning = False
        self.pool = QThreadPool()
        self.savename = ''
        self.savefolder = str(Path.home())
        self.setup_flag=False

        super().__init__()


    def connect(self, d0=None) -> None:

        self.api.connect()

        time.sleep(1)
        self.setup(d0)
        
        return super().connect()

    def setup(self,d0=None):

        self.api.hidden['BrPlsMode'].value = True
        self.api.hidden['OpMode'].value = 'Operate'
        # self.api.hidden['RefArm'].value = 'On'

        self.api.cur_exp['VideoBW'].value = 20
        if d0 is None:
            self.calc_d0()
        else:
            self.d0 = d0
        self.api.hidden['Detection'].value = 'Signal'
        self.setup_flag = True
        pass

    def acquire_dataset(self):
        if self.bg_data is None:
            if not self.isrunning():
                if (self.savename is not None) and (self.savename != ''):
                    self.api.xepr_save(os.path.join(self.savefolder,self.savename))
                data = self.api.acquire_dataset(self.cur_exp)
            else:
                data = self.api.acquire_scan(self.cur_exp)
        else:

            data = create_dataset_from_sequence(self.bg_data, self.cur_exp)
        return super().acquire_dataset(data)


    def _launch_complex_thread(self,sequence,axID=1,tune=True):
    
        uProgTable = build_unique_progtable(sequence)
        reduced_seq = sequence.copy()
        reduced_seq.progTable = transpose_list_of_dicts([transpose_dict_of_list(sequence.progTable)[0]])

        uProgTable_py = uProgTable[axID]
        axis = uProgTable_py['axis']
        reduced_seq.averages.value = 1
        py_ax_dim = uProgTable_py['axis']['dim']

        
        self.bg_data = np.zeros((uProgTable[0]['axis']['dim'],py_ax_dim),dtype=np.complex128)
        
        print("Initial PulseSpel Launch")
        self.launch(reduced_seq,savename='test',tune=tune, update_pulsespel=True, start=False,reset_bg_data=False,reset_cur_exp=False)

        self.terminate()

        variables = uProgTable_py['variables']
        # print("Creating Thread")
        # # thread = threading.Thread(target=step_parameters,args=[self,reduced_seq,py_ax_dim,variables])
        # # self.bg_thread = self.pool.submit(step_parameters, self,reduced_seq,py_ax_dim,variables)
        # self.bg_thread = Worker(step_parameters, self,reduced_seq,py_ax_dim,variables)
        # self.pool.start(self.bg_thread)
        # print("Started Thread")

        step_parameters(self,reduced_seq,py_ax_dim,variables)
        # thread.start()

        pass
        
    def launch(self, sequence: Sequence, savename: str, start=True, tune=True,
            MPFU_overwrite=None,update_pulsespel=True, reset_bg_data = True,
            reset_cur_exp=True,**kwargs):
        
        sequence.shift_detfreq_to_zero()

        if self.isrunning():
            self.terminate(now=True)
            time.sleep(4)
                    
        if reset_bg_data:
            self.bg_data = None

        if reset_cur_exp:
            timestamp = datetime.datetime.now().strftime(r'%Y%m%d_%H%M_')
            self.savename = timestamp+savename
            self.cur_exp = sequence
                    
        # First check if the sequence is pulsespel compatible
    
        if not test_if_MPFU_compatability(sequence):
            print("Launching complex sequence")
            self._launch_complex_thread(sequence,1,tune)
            return None
                
        
        if update_pulsespel:
            
            def_text, exp_text, pulse_shapes = write_pulsespel_file(sequence,self.d0,AWG=True)
            
            verbMsgParam = self.api.cur_exp.getParam('*ftEPR.PlsSPELVerbMsg')
            plsSPELCmdParam = self.api.cur_exp.getParam('*ftEPR.PlsSPELCmd')
            self.api.XeprCmds.aqPgSelectBuf(1)
            self.api.XeprCmds.aqPgSelectBuf(2)
            self.api.XeprCmds.aqPgSelectBuf(3)

            # Merge the pulse shapes into one string
            pulse_shapes = '\n'.join(pulse_shapes)

            self.api.cur_exp.getParam('*ftEpr.PlsSPELShpTxt').value = pulse_shapes

            plsSPELCmdParam.value=9
            time.sleep(2)
            self.api.XeprCmds.aqPgSelectBuf(1)

            self.api.cur_exp.getParam('*ftEpr.PlsSPELGlbTxt').value = def_text
            self.api.XeprCmds.aqPgShowDef()

            plsSPELCmdParam.value=3
            time.sleep(2)
            # while not "The variable values are set up" in verbMsgParam.value:
            #     time.sleep(0.1)

            self.api.XeprCmds.aqPgSelectBuf(2)

            self.api.cur_exp.getParam('*ftEpr.PlsSPELPrgTxt').value = exp_text
            plsSPELCmdParam.value=7
            time.sleep(2)
            # while not "Second pass ended" in verbMsgParam.value:
            #     time.sleep(0.1)


                
        self.api.set_field(sequence.B.value)
        self.api.set_freq(sequence.LO.value)
        
        if 'B' in sequence.progTable['Variable']:
            idx = sequence.progTable['Variable'].index('B')
            B_axis = sequence.progTable['axis'][idx]
            self.api.set_sweep_width(B_axis.max()-B_axis.min())
        

        self.api.set_ReplaceMode(False)
        self.api.set_Acquisition_mode(1)
        self.api.set_PhaseCycle(True)
        pg = sequence.pulses[-1].tp.value
        pg = round_step(pg/2,self.bridge_config['Pulse dt'])
        d0 = self.d0-pg

        if d0 <0:
            d0=0
        d0 = round_step(d0,self.bridge_config['Pulse dt'])
        self.api.set_PulseSpel_var('d0', d0)
        self.api.set_PulseSpel_var('pg', pg*2)
        self.api.set_PulseSpel_experiment('auto')
        self.api.set_PulseSpel_phase_cycling('auto')

        if start:
            self.api.run_exp()
        pass

    
    def tune_rectpulse(self,*,tp, LO, B, reptime, shots=400):
        """Generates a rectangular pi and pi/2 pulse of the given length at 
        the given field position. This value is stored in the pulse cache. 

        Parameters
        ----------
        tp : float
            Pulse length in ns
        LO : float
            Central frequency of this pulse in GHz
        B : float
            Magnetic B0 field position in Gauss
        reptime: float
            Shot repetion time in us.
        shots: int
            The number of shots

        Returns
        -------
        p90: RectPulse
            A tuned rectangular pi/2 pulse of length tp
        p180: RectPulse
            A tuned rectangular pi pulse of length tp
        """

        amp_tune =HahnEchoSequence(
            B=B, LO=LO, reptime=reptime, averages=1, shots=shots
        )

        scale = Parameter("scale",0,dim=51,step=0.02)
        amp_tune.pulses[0].tp.value = tp
        amp_tune.pulses[0].scale = scale
        amp_tune.pulses[1].tp.value = tp * 2
        amp_tune.pulses[1].scale = scale

        amp_tune.evolution([scale])

        vg = self.api.get_video_gain()
        overflow_flag= True
        while overflow_flag:
            
            self.launch(amp_tune, "autoDEER_amptune", IFgain=1)

            while self.isrunning():
                time.sleep(10)
            dataset = self.acquire_dataset()
            dataset = dataset.epr.correctphase

            data = np.abs(dataset.data)
            scale_amp = np.around(dataset.pulse0_scale[data.argmax()].data,2)
            if scale_amp > 0.95:
                raise RuntimeError("Not enough power avaliable.")
            
            if scale_amp == 0:
                warnings.warn("Pulse tuned with a scale of zero!")
            p90 = amp_tune.pulses[0].copy(
                scale=scale_amp, LO=amp_tune.LO)
            
            p180 = amp_tune.pulses[1].copy(
                scale=scale_amp, LO=amp_tune.LO)
            

            scale = Parameter("scale",scale_amp,dim=4,step=0.01)
            amp_tune.pulses[0].scale = scale
            amp_tune.pulses[1].scale = scale
            amp_tune.evolution([scale])

            self.launch(amp_tune, "autoDEER_amptune", IFgain=1)

            time.sleep(3)
            self.terminate()
            self.api.cur_exp['ftEPR.StartPlsPrg'].value = True
            self.api.hidden['specJet.NoOfPoints'].value = 512
            self.api.hidden['specJet.NoOfAverages'].value = 20
            self.api.hidden['specJet.NOnBoardAvgs'].value = 20
            if not self.api.hidden['specJet.AverageStart'].value:
                self.api.hidden['specJet.AverageStart'].value = True
            
            limit = np.abs(self.api.hidden['specjet.DataRange'][0])
            data = get_specjet_data(self)
            
            if (np.abs(data.real).max() > 0.7*limit) or (np.abs(data.imag).max() > 0.7*limit):
                self.api.set_video_gain(vg - 9)
                overflow_flag= True
                print('overflow')
            elif (np.abs(data.real).max() < 0.3*limit) and (np.abs(data.imag).max() < 0.3*limit):
                self.api.set_video_gain(vg + 6)
                overflow_flag= True
                print('underflow')
            else:
                overflow_flag = False

        return p90, p180
    
    def tune_pulse(self, pulse, mode, LO, B , reptime, shots=400):
        """Tunes a single pulse a range of methods.

        Parameters
        ----------
        pulse : Pulse
            The Pulse object in need of tuning.
        mode : str
            The method to be used.
        LO : float
            The local oscilator frequency in GHz
        B : float
            Magnetic B0 field position in Gauss
        reptime : us
            Shot repetion time in us.
        shots: int
            The number of shots

        Returns
        -------
        Tunned Pulse: Pulse
            The returned pulse object that is now tunned.
        """
        # Check pulse is a pulse
        if type(pulse) == Delay:
            pass
        if type(pulse) == Detection:
            pass
        
        # Get absolute central frequency
        if hasattr(pulse,"freq"):
            c_frq = pulse.freq.value + LO
        elif hasattr(pulse, "init_freq") & hasattr(pulse, "BW"):
            c_frq = pulse.init_freq.value + 0.5*pulse.BW.value + LO
        elif hasattr(pulse, "final_freq") & hasattr(pulse, "BW"):
            c_frq = pulse.final_freq.value - 0.5*pulse.BW.value + LO
        elif hasattr(pulse, "init_freq") & hasattr(pulse, "final_freq"):
            c_frq = 0.5*(pulse.final_freq.value + pulse.final_freq.value) + LO

        pi2_pulse, pi_pulse = self.tune_rectpulse(tp=12, B=B, LO=c_frq, reptime=reptime)
        
        if mode == "amp_hahn":
            amp_tune =HahnEchoSequence(
                B=B, LO=LO, 
                reptime=reptime, averages=1, shots=shots,
                pi2_pulse = pulse, pi_pulse=pi_pulse
            )

            scale = Parameter('scale',0,unit=None,step=0.02, dim=51, description='The amplitude of the pulse 0-1')
            amp_tune.pulses[0].scale = scale

            amp_tune.evolution([scale])

            self.launch(amp_tune, "autoDEER_amptune", IFgain=1)

            while self.isrunning():
                time.sleep(10)
            dataset = self.acquire_dataset()
            new_amp = np.around(dataset.pulse0_scale[dataset.data.argmax()].data,2)
            pulse.scale = Parameter('scale',new_amp,unit=None,description='The amplitude of the pulse 0-1')
            return pulse

        elif mode == "amp_nut":

            nut_tune = Sequence(
                name="nut_tune", B=(B/LO*c_frq), LO=LO, reptime=reptime,
                averages=1,shots=shots
            )
            nut_tune.addPulse(pulse.copy(
                t=0, pcyc={"phases":[0],"dets":[1]}, scale=0))
            nut_tune.addPulse(
                pi2_pulse.copy(t=2e3,
                               pcyc={"phases":[0, np.pi],"dets":[1, -1]},
                               freq=c_frq-LO))
            nut_tune.addPulse(
                pi_pulse.copy(t=2.5e3, pcyc={"phases":[0],"dets":[1]},
                              freq=c_frq-LO))
            nut_tune.addPulse(Detection(t=3e3, tp=512, freq=c_frq-LO))

            scale = Parameter('scale',0,unit=None,step=0.02, dim=51, description='The amplitude of the pulse 0-1')
            nut_tune.pulses[0].scale = scale
            nut_tune.evolution([scale])

            self.launch(nut_tune, "autoDEER_amptune", IFgain=1)

            while self.isrunning():
                time.sleep(10)
            dataset = self.acquire_dataset()
            dataset = dataset.epr.correctphase
            data = dataset.data
            axis = dataset.scale
            # data = correctphase(dataset.data)
            if data[0] < 0:
                data *= -1

            if np.isclose(pulse.flipangle.value, np.pi):
                new_amp = np.around(axis[data.argmin()].data,2)
            elif np.isclose(pulse.flipangle.value, np.pi/2):
                sign_changes = np.diff(np.sign(np.real(data)))
                new_amp = np.around(axis[np.nonzero(sign_changes)[0][0]].data,2)
            else:
                raise RuntimeError("Target pulse can only have a flip angle of either: ",
                                "pi or pi/2.")

            pulse.scale = Parameter('scale',new_amp,unit=None,description='The amplitude of the pulse 0-1')
        
            return pulse

    def phasetune_pulse(self, pulse):
        # Bruker SpecJet-I has a phase issue. The pulse is observed and optimised through the transmision mode.

        if isinstance(pulse, Detection):
            raise RuntimeError("Detection pulses cannot be phase tuned.")
        
        current_B = self.api.get_field()
        current_LO = self.api.get_counterfreq()

        test_seq = Sequence(name='test_seq',B=current_B,LO=current_LO,reptime=250,averages=1,shots=20)
        test_seq.addPulse(pulse)
        test_seq.addPulse(Detection(tp=16,t=-self.d0+100))



        
    def isrunning(self) -> bool:
        if self.tuning:
            return True
        if self.bg_thread is None:
            return self.api.is_exp_running()
        else:
            return self.bg_thread.running()

    def terminate(self,now=False) -> None:
        self.tuning = False
        if self.bg_thread is None:
            if now:
                return self.api.abort_exp()
            else:
                return self.api.stop_exp()
        else:
            attempt = self.bg_thread.cancel()
            if not attempt:
                raise RuntimeError("Thread failed to be canceled!")
            
    def calc_d0(self):
        """
        This creates an initial guess for d0. A better estimate can only be found after the field sweep. 
        """
        hw_log.info('Calcuating d0')
        hw_log.debug('Setting Detection = TM')
        self.api.hidden['Detection'].value = 'TM'
        B = self.api.get_field()
        LO = self.api.get_counterfreq()

        self.api.set_attenuator('+<x>',100)

        d0=0
        self.d0=d0

        seq = Sequence(name='single_pulse',B=B,LO=LO,reptime=3e3,averages=1,shots=20)
        det_tp = Parameter('tp',value=16,dim=4,step=0)
        seq.addPulse(RectPulse(tp=det_tp,t=0,flipangle=np.pi,scale=1))
        seq.addPulse(Detection(tp=16,t=d0))

        seq.evolution([det_tp])
        self.launch(seq,savename='test',tune=False)
        self.terminate(now=True)
        time.sleep(3)

        self.api.cur_exp['ftEPR.StartPlsPrg'].value = True
        self.api.hidden['specJet.NoOfAverages'].value = 20
        self.api.hidden['specJet.NOnBoardAvgs'].value = 20
                
        if not self.api.hidden['specJet.AverageStart'].value:
            self.api.hidden['specJet.AverageStart'].value = True

        self.api.hidden['specJet.NoOfPoints'].value = 1024
        time.sleep(3)
        
        
        optimal = False
        while not optimal:
            max_value = np.abs(get_specjet_data(self)).max()
            y_max = np.abs(self.api.hidden['specjet.DataRange'][0])
            vg  =self.api.get_video_gain()
            if max_value > 0.7* y_max:
                self.api.set_video_gain(vg - 3)
                time.sleep(0.5)
            elif max_value < 0.3* y_max:
                self.api.set_video_gain(vg + 3)
                time.sleep(0.5)
            else:
                optimal=True
        
        specjet_data = np.abs(get_specjet_data(self))

        calc_d0 = d0  + self.api.hidden['specJet.Abs1Data'][specjet_data.argmax()]

        d0 = calc_d0 - 256

        seq = Sequence(name='single_pulse',B=B,LO=LO,reptime=3e3,averages=1,shots=20)
        det_tp = Parameter('tp',value=16,dim=4,step=0)
        seq.addPulse(RectPulse(tp=det_tp,t=0,flipangle=np.pi,scale=1))
        seq.addPulse(Detection(tp=16,t=d0))

        seq.evolution([det_tp])
        self.launch(seq,savename='test',tune=False)
        self.terminate()

        self.api.cur_exp['ftEPR.StartPlsPrg'].value = True
                
        if not self.api.hidden['specJet.AverageStart'].value:
            self.api.hidden['specJet.AverageStart'].value = True

        self.api.hidden['specJet.NoOfPoints'].value = 512
        time.sleep(3)
        specjet_data = np.abs(get_specjet_data(self))
        
        calc_d0 = d0 + self.api.hidden['specJet.Abs1Data'][specjet_data.argmax()]

        self.d0 = calc_d0 + 64 # 64ns added to compensate for hahn echo center in field sweep
        
        hw_log.info(f"d0 set to {self.d0}")
        self.api.hidden['Detection'].value = 'Signal'

    def calc_d0_from_Hahn_Echo(self, B=None, LO=None):
        
        B = self.api.get_field()
        LO = self.api.get_counterfreq()

        if B is not None:
            self.api.set_field(B)
        if LO is not None:
            self.api.set_freq(LO)
        
        d0 = self.d0
        # self.api.set_PulseSpel_var('d0',d0)
        self.api.run_exp()
        self.api.abort_exp()
        while self.api.is_exp_running():
            time.sleep(1)

        self.api.cur_exp['ftEPR.StartPlsPrg'].value = True
        self.api.hidden['specJet.NoOfAverages'].value = 20
        self.api.hidden['specJet.NOnBoardAvgs'].value = 20
                
        if not self.api.hidden['specJet.AverageStart'].value:
            self.api.hidden['specJet.AverageStart'].value = True

        self.api.hidden['specJet.NoOfPoints'].value = 512

        optimal = False
        while not optimal:
            max_value = np.abs(get_specjet_data(self)).max()
            y_max = np.abs(self.api.hidden['specjet.DataRange'][0])
            vg  =self.api.get_video_gain()
            if max_value > 0.7* y_max:
                self.api.set_video_gain(vg - 3)
                time.sleep(0.5)
            elif max_value < 0.3* y_max:
                self.api.set_video_gain(vg + 3)
                time.sleep(0.5)
            else:
                optimal=True
        
        specjet_data = np.abs(get_specjet_data(self))
        calc_d0 = d0 - 64 + self.api.hidden['specJet.Abs1Data'][specjet_data.argmax()]

        self.d0 = calc_d0 
        hw_log.info(f"d0 set to {self.d0}")
        return self.d0

# =============================================================================


def get_specjet_data(interface):
    n_points  = interface.api.hidden['specJet.Data'].aqGetParDimSize(0)
    array = np.zeros(n_points,dtype=np.complex128)
    for i in range(n_points):
        array[i] = interface.api.hidden['specJet.Data'][i] + 1j* interface.api.hidden['specJet.Data1'][i]

    return array

def test_if_MPFU_compatability(seq):
    table = seq.progTable
    if 'LO' in table['Variable']:
        return False
    elif np.unique(table['axID']).size > 2:
        return False 
    else:
        return True