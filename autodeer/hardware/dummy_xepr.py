# -*- coding: utf-8 -*-
"""
This is a set of scripts that pretend to be an xepr spectromter, this allows for the further development of 
this code without having to be sat in front of an actual spectrometer:

Manufactuer: Hugo Karas (me)
Model: DUMMY
Features: The most relibable Bruker EPR spectromter in the world (possibly?)
"""

import re
import numpy as np
import deerlab as dl
import time
import os
import random as rand
from .XeprAPI_link import XeprAPILink,dataset

import logging

hw_log = logging.getLogger('hardware.dummy')

rand.seed(1212)

hardware_meta = {# This dictionary should be moved into a config file eventually
    "Type":             "Complete Spectrometer",
    "Manufactuer":      "Hugo Karas",
    "Model":            "Dummy",
    "Local name":       "Dummy",
    "Acq. Resolution":  2,
    "Pulse Resolution": 2,
    "AWG":              False,
    "Min Freq":         33,
    "Max Freq":         35,
 } 

class dummy_api(XeprAPILink):

    def __init__(self) -> None:
        self._tmp_dir = None
        self.meta = hardware_meta # This will become more neuanced eventually. Currentlty this is sufficent.
        pass

    def find_Xepr(self):
        self.Xepr = dummy_xepr()
        pass

    def is_exp_running(self):
        # Since all experiments are instantaneous this has to be always false
        return False

    def acquire_scan(self):
        return self.acquire_dataset()
    
    def acquire_scan_at(self, scan_num: int):
        return self.acquire_dataset()
    
    def acquire_scan_2d(self):
        return self.acquire_dataset()

    def compile_PulseSpel_def(self):
        pass

    def compile_PulseSpel_prg(self):
        pass

    def set_PulseSpel_var(self, variable: str, value: int):
        # This should be expanded into a dictionary that can then effect the simulation. In an ideal world
        pass
    def run_exp(self):
        hw_log.info('Experiment started')
        pass

    def stop_exp(self):
        hw_log.info('Experiment stopped')
        pass
    
    def abort_exp(self):
        hw_log.info('Experiment aborted')
        pass
    

class dummy_cur_exp:
    
    def __init__(self):
        # if experiment == "DEER_quick":
        #     print("Experiment set to: DEER_quick")
        #     self.NbScansToDo = dummy_param(0,10000,int,112)
        #     self.SweepsPExp = dummy_param(0,10000,int,112)                      # This is the same as scans in a 1D experiment
        #     self.ShotsPLoop = dummy_param(0,10000,int,4)                        # This is the same as shots per point
        #     self.XSpecRes   = dummy_param(0,10000,int,256)                      # This is the length of the X direction
        #     self.NbScansDone = dummy_param(0,self.NbScansToDo,int,0)            # Intially zero scans have been done 
        #     self.FTAcqModeSlct = dummy_param(None,None,str,'Run from PulseSPEL')# Setting teh acquistion mode
        #     self.PlsSPELPrgPaF = dummy_param(None,None,str,'rand_path')         # Pulse Spel Experiment Path
        #     self.PlsSPELGlbPaF = dummy_param(None,None,str,'rand_path')         # Pulse Spel Definitions Path
        #     self.PlsSPELLISTSlct = dummy_param(None,None,str,'Phase Cycle')     # Selecting the phase cycling
        #     self.PlsSPELEXPSlct = dummy_param(None,None,str,'Experiment')       # Selecting the Experiment
        #     self.ReplaceMode = dummy_param(None,None,bool,False)                # Setting the replace mode state

        # elif experiment == "DEER_std":#TODO
        #     print("Experiment set to: DEER_std")
        #     self.NbScansToDo = dummy_param(0,10000,int,2000)
        #     self.SweepsPExp = dummy_param(0,10000,int,112)
        # elif experiment == "2D_DEC":#TODO
        #     print("Experiment set to 2D Dec")   

        self.ShotRepTime = dummy_param(0,100000,int,6000)
        self.NbScansToDo = dummy_param(0,10000,int,112)
        self.SweepsPExp = dummy_param(0,10000,int,112)                      # This is the same as scans in a 1D experiment
        self.ShotsPLoop = dummy_param(0,10000,int,4)                        # This is the same as shots per point
        self.XSpecRes   = dummy_param(0,10000,int,256)                      # This is the length of the X direction
        self.NbScansDone = dummy_param(0,self.NbScansToDo,int,0)            # Intially zero scans have been done 
        self.FTAcqModeSlct = dummy_param(None,None,str,'Run from PulseSPEL')# Setting teh acquistion mode
        self.PlsSPELPrgPaF = dummy_param(None,None,str,'rand_path')         # Pulse Spel Experiment Path
        self.PlsSPELGlbPaF = dummy_param(None,None,str,'rand_path')         # Pulse Spel Definitions Path
        self.PlsSPELLISTSlct = dummy_param(None,None,str,'Phase Cycle')     # Selecting the phase cycling
        self.PlsSPELEXPSlct = dummy_param(None,None,str,'Experiment')       # Selecting the Experiment
        self.ReplaceMode = dummy_param(None,None,bool,False)                # Setting the replace mode state
        self.CenterField = dummy_param(0,15000,int)
        self.FrequencyMon = dummy_param(33,35,float)
        self.SweepWidth = dummy_param(0,10000,int,300)
    

    def __getitem__(self, name):
        """This method returns a given parameter class for the name provided

        :param name: The name of parameter needed
        :type name: string
        :return: Instance of dummy_param class corresponding to the searched name
        :rtype: dummy_param class
        """
        if '.' in name:
            # This removes the category allowing the paramter class to be found
            name = name.split('.')[1]
        return(getattr(self,name))

    def getParam(self,name:str):
        """Extracts a paremeter from the current exeperiment

        :param name: The name of the parameter to be extracted
        :type name: string
        :return: An instance of the parameter class for the given paramter name
        :rtype: parameter Class
        """
        return getattr(self,name)

    def aqExpRun(self):
        #Dummy function, for running the experiment
        print('Dummy Experiment running')
        pass


class dummy_param:

    def __init__(self,min,max,type=int,par=None):
        self.aqGetParMaxValue = max
        self.aqGetParMinValue = min
        self.val = par 
        if (par == None) & (type == int): # If no value is given a random one is picked
            self.val = rand.randint(self.aqGetParMinValue,self.aqGetParMaxValue)

    @property
    def value(self):
        return self.val

    @value.setter
    def value(self, par):
        self.val = par

class dummy_xepr:

    def __init__(self) -> None:
        self.cur_exp = dummy_cur_exp()
        self.hidden = dummy_hidden()
        pass

    def XeprExperiment(self,hidden=None):
        if hidden == None:
            # return dummy_cur_exp("DEER_quick")
            return self.cur_exp
        elif hidden == "AcqHidden":
            return self.hidden
    
    def XeprDataset(self): # Returns the current dataset
        
        cur_exp = self.cur_exp
        hidden = self.hidden
        # Generates simulated DEER data
        def generateDEER(num_points,max_time):
            """
            Modified from an example in DEERlab.

            Copyright 2019-2021, Luis Fábregas Ibáñez, Stefan Stoll, and others.
            """
            t = np.linspace(0,max_time/1000,num_points)        # time axis, µs
            r = np.linspace(2,5,200)           # distance axis, nm
            param = [3, 0.1, 0.2, 3.5, 0.1, 0.65, 3.8, 0.05, 0.15] # parameters for three-Gaussian model
            P = dl.dd_gauss3(r,param)          # model distance distribution
            lam = 0.3                          # modulation depth
            B = dl.bg_hom3d(t,300,lam)         # background decay
            K = dl.dipolarkernel(t,r,mod=lam,bg=B)    # kernel matrix
            Vexp = K@P + dl.whitegaussnoise(t,0.01,seed=0)  # DEER signal with added noise
            return t, Vexp

        def generateCarrPurcell(num_points,max_time):
            t = np.linspace(0,4000,256)
            V = np.exp(-((t*3e-3)/(4.2))**4)
            noise = np.random.normal(0, .05, V.shape) *np.array([1+1j]) 
            Vexp = V + noise
            return t, Vexp

        def generate2D(numpoints,max_time):
            t = np.linspace(0,max_time,numpoints)
            X, Y = np.meshgrid(t, t)
            def fun_2d(x,y):
                return (np.exp(-((x*2.8e-3)/(4.2))**4) + np.exp(-((y*2.8e-3)/(4.2))**4))
            Z = fun_2d(X,Y)
            noise = np.random.normal(0, .05, Z.shape) *np.array([1+1j])
            return [t,t],Z+noise

        def generate_hahn_echo_transient(numpoints,phase,phase_offset=0.2):
            def gaussian(t,a,b,c,phase):
                return a * np.exp(-(t-b)**2/c) * np.exp(-1j*(phase+phase_offset))
            t = np.arange(0,numpoints,1) 
            noise = np.random.normal(loc=0.0,scale=0.1,size=t.shape) + 1j * np.random.normal(loc=0.0,scale=0.03,size=t.shape)
            V = gaussian(t,1,numpoints/2,numpoints*2,phase) + noise
            return t,V

        def generate_phase_experiment(numpoints,pg,phase,phase_offset=0.2):
            vals = np.zeros(numpoints)
            vals = vals + 1j* vals
            t = range(0,numpoints)
            for i in t:
                ta,V = generate_hahn_echo_transient(pg,phase,phase_offset)
                vals[i,] = np.trapz(V,x=ta)
            return t, vals

        
        if cur_exp.PlsSPELEXPSlct.value == 'DEER':
            t,V = generateDEER(250,4000)
            return dummy_dataset(t,V)
        elif cur_exp.PlsSPELEXPSlct.value =='2D Dec. 64':
            t,V = generate2D(64,4000)
            return dummy_dataset(t,V)
        elif cur_exp.PlsSPELEXPSlct.value =='Carr Purcell exp':
            t,V = generateCarrPurcell(250,4000)
            return dummy_dataset(t,V)
        elif cur_exp.PlsSPELEXPSlct.value =='tau 2 scan':
            t,V = generateCarrPurcell(250,4000)
            return dummy_dataset(t,V)
        elif cur_exp.PlsSPELEXPSlct.value == 'Hahn Echo':
            def phases(x,range):
                m = 6/range *np.pi
                c = -3*np.pi
                return m*x + c
            if cur_exp.PlsSPELLISTSlct.value == 'MainPhase':
                phase_setting = hidden['SignalPhase'].value
                phase = phases(phase_setting,hidden['SignalPhase'].aqGetParMaxValue)
            elif cur_exp.PlsSPELLISTSlct.value in ['BrXPhase','BrYPhase','MinBrXPhase','MinBrYPhase']:
                phase_setting = hidden[cur_exp.PlsSPELLISTSlct.value].value
                phase = phases(phase_setting,hidden[cur_exp.PlsSPELLISTSlct.value].aqGetParMaxValue)
            else:
                print("Sorry, This can not be simulated yet, returning Main Phase")
                phase_setting = hidden['SignalPhase']
                phase = phases(phase_setting,hidden['SignalPhase'].aqGetParMaxValue())
            
            t,V = generate_phase_experiment(10,512,phase)
            return dummy_dataset(t,V)
        else:
            print("Sorry, This can not be simulated yet, returning white noise")
            t = np.linspace(0,2500,250)
            noise = np.random.normal(0, .1, t.shape) *np.array([1+1j])
            return dummy_dataset(t,noise)        
        
    class XeprCmds:

        def aqPgShowPrg():
            # Opens the Pulse Spel Panel, obvs does nothing
            pass
        def aqPgCompile():
            # Compiles Pulse Spel, obvs does nothing
            pass

        def aqExpSelect(exp):
            # Would normally select the experiment, dummy does nothing
            pass

        def aqPgLoad(filepath):
            # Loads the pulse spel program, does nothing
            pass
        
        def aqPgDefLoad(filepath):
            # Loas the pulse spel definitions file, does nothing
            pass

        def aqPgDefSaveAs(path):
            pass

        def aqPgSaveAs(path):
            pass
        
        def vpSave(viewpoint:str,title:str,path:str):
            pass

class dummy_hidden:
    # This is a dummy version of the hidden class containing the extra current experiment features

    def __init__(self) -> None:
        self.SignalPhase = dummy_param(0,4095,type=int,par=rand.randint(0,4096))
        self.BrXPhase = dummy_param(0,100,type=float,par=rand.random())
        self.BrYPhase = dummy_param(0,100,type=float,par=rand.random())
        self.MinBrXPhase = dummy_param(0,100,type=float,par=rand.random())
        self.MinBrYPhase = dummy_param(0,100,type=float,par=rand.random())

        self.BrXAmp = dummy_param(0,100,type=float,par=rand.random())
        self.BrYAmp = dummy_param(0,100,type=float,par=rand.random())
        self.MinBrXAmp = dummy_param(0,100,type=float,par=rand.random())
        self.MinBrYAmp = dummy_param(0,100,type=float,par=rand.random())
       
        pass


    def __getitem__(self, name):
        """This method returns a given parameter class for the name provided

        :param name: The name of parameter needed
        :type name: string
        :return: Instance of dummy_param class corresponding to the searched name
        :rtype: dummy_param class
        """
        if '.' in name:
            # This removes the category allowing the paramter class to be found
            name = name.split('.')[1]
        return(getattr(self,name))

    def getParam(self,name:str):
        """Extracts a paremeter from the current exeperiment

        :param name: The name of the parameter to be extracted
        :type name: string
        :return: An instance of the parameter class for the given paramter name
        :rtype: parameter Class
        """
        return getattr(self,name)

class dummy_dataset:

    def __init__(self,t,V) -> None:
        self.O = V
        self.size = np.shape(V)
        if len(self.size) == 1:
            self.X = t
        elif len(self.size) == 2:
            self.X = t[0]
            self.Y = t[1]
    
