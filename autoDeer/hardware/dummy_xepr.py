""""
This is a set of scripts that pretend to be an xepr spectromter, this allows for the further development of 
this code without having to be sat in front of an actual spectrometer:

Manufactuer: Hugo Karas (me)
Model: DUMMY
Features: The most relibable Bruker EPR spectromter in the world (possibly?)
"""

import numpy as np
import time
import os
import random as rand

rand.seed(1212)

Xepr = None
cur_exp = None
hidden = None

def get_Xepr_global():
    """This function gets the Xepr class

    :raises RuntimeError: This error is raised if the spectometer has not been started
    :return: The current Xepr instance 
    :rtype: [type]
    """
    if Xepr != None:
        return Xepr
    else:
        raise RuntimeError("Please intialize the dummy spectrometer")


def find_cur_exp(experiment:str):

    if experiment == "DEER_quick":
        print("Experiment set to: DEER_quick")
    
    elif experiment == "DEER_std":
        print("Experiment set to: DEER_std")
    
    elif experiment == "2D_DEC":
        print("Experiment set to 2D Dec")

    # These are just general random settings

class dummy_cur_exp:
    
    def __init__(self,experiment:str):
        if experiment == "DEER_quick":
            print("Experiment set to: DEER_quick")
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

        elif experiment == "DEER_std":#TODO
            print("Experiment set to: DEER_std")
            self.NbScansToDo = dummy_param(0,10000,int,2000)
            self.SweepsPExp = dummy_param(0,10000,int,112)
        elif experiment == "2D_DEC":#TODO
            print("Experiment set to 2D Dec")   

        self.ShotRepTime = dummy_param(0,100000,int,6000)
        
    

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

    def aqExpRun():
        #Dummy function, for running the experiment
        print('Dummy Experiment running')
        pass


class dummy_param:

    def __init__(self,min,max,type=int,value=None):
        self.aqGetParMaxValue = max
        self.aqGetParMinValue = min
        if value == None & type == int: # If no value is given a random one is picked
            self.value = rand.randint(self.aqGetParMinValue,self.aqGetParMaxValue)

    @property
    def value(self):
        return self.value

    @value.setter
    def value(self, val):
        self.value = val

class dummy_xepr:

    def XeprExperiment(hidden=None):
        if hidden == None:
            return dummy_cur_exp("DEER_quick")
        elif hidden == "AcqHidden":
            return dummy_hidden
    
    
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

class dummy_hidden:
    # This is a dummy version of the hidden class containing the extra current experiment features

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