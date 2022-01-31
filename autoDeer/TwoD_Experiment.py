# This file contains the class for running 2D Decoherence experiments
# all code assoicated with processing these experiments are included in
# this file.

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
from deerlab import deerload
import os

class TwoD_Experiment:
    """"
    TwoD_Experiment: This is a class for loading and processing 2D Decoherence experiments
    """

    def __init__(self) -> None:
        self.snr_target = 20 #normal units not db
        self.time_target = 2 #hrs
        self.trace_length = 512 #this is wrong, and is likely to be a function of the optimal time. 
        self.noise_frac = 0.75


    def load(self,file):
        """
        Using deerlab the file is loaded. If the file is a bruker .DSC file then extra paremeters are loaded. 
        If it is a .csv or .txt file then any required metadata will need to be inputted by the user
        Parameters:
        self:
        file (string): This is the full file path for the datafile
        """
        file_name,file_ext = os.path.splitext(file)
        
        if file_ext == '.DSC':
            self.time,self.data,self.params = deerload(file,full_output=True)
            self.data = self.data.T
            self.scans = int(self.params.get('DSL').get('recorder').get('NbScansDone'))
            self.shots = int(self.params.get('DSL')['ftEpr']['ShotsPLoop'])
            self.shrt = float(self.params.get('DSL')['ftEpr']['ShotRepTime'].split()[0]) # in us
        else:
            self.time,self.data= deerload(file,full_output=False)
            self.params = None
            
    def import_data(self,time,data,scans:int,shots:int,shrt:int):
        """
        The loads all the import infomation directly from python variables and arrays
        Parameters:
        time - [timex,timey]
        data - numpy 2d array
        scans(int)
        shots(int)
        shrt(int)
        """
        self.time = time
        self.data = data
        self.scan = scans
        self.shots = shots
        self.shrt = shrt

    def import_dataset(self,dataset):
        self.time = dataset.time
        self.data = dataset.data
        self.scan = dataset.scans_done
        self.shots = dataset.shot_p_point
        self.shrt = dataset.shrt

    def create_bahrenberg_plots(self):
        """
        Returns a matplotlib figure object for the standard bahrenberg figure.
        """
        
        time = self.time
        signal_real = np.real(self.data)
        cmap = cm.get_cmap(name='viridis', lut=None)
        cmap2 = cm.get_cmap(name='Greys', lut=None)
        fig = plt.figure(figsize=(12, 6),dpi=150)
        axs = fig.subplots(1,2)
        
        signal_real = signal_real - signal_real.min()
        signal_real = signal_real/ signal_real.max()
        scale = [min(time[0]),max(time[0]),min(time[1]),max(time[1])]

        tau2_max_locs = signal_real.argmax(axis=1)
        tau1_max_locs = signal_real.argmax(axis=0)

        axs[0].contour(signal_real,extent=scale,levels=10,cmap=cmap2)
        axs[0].pcolormesh(time[0],time[1],signal_real,shading='auto',alpha=0.95,cmap=cmap)
        axs[0].plot([min(time[0]),max(time[0])],[min(time[1]),max(time[1])],color='k',linestyle='--')
        axs[0].plot(time[0][tau2_max_locs],time[1],color='r')
        axs[0].plot(time[0],time[1][tau2_max_locs],color='r')
        axs[0].set_xlabel(r'$\tau_1 / (us)$')
        axs[0].set_ylabel(r'$\tau_2 / (us)$')
        axs[0].set_aspect('equal')


        # Second Tau2 Normalised Plot
        signal_real_tau2_norm = signal_real / signal_real.max(axis=1)[:,None]
        im = axs[1].pcolormesh(time[0],time[1],signal_real_tau2_norm,shading='auto',alpha=0.95,cmap=cmap)

        axs[1].contour(time[0],time[1],signal_real_tau2_norm,levels=5,cmap=cmap2)
        axs[1].plot([min(time[0]),max(time[0])],[min(time[1]),max(time[1])],color='k',linestyle='--')
        axs[1].plot(time[0][tau2_max_locs],time[1],color='r')
        axs[1].set_xlabel(r'$\tau_1 / (us)$')
        axs[1].set_ylabel(r'$\tau_2 / (us)$')
        axs[1].set_aspect('equal')


        cb_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])

        fig.colorbar(im,cax=cb_ax)

        return fig
    
    def _stability_check(self):
        """Private method that caclulates if a value in the 2D decoherence experiment is above the snr_threshold. 
        The complete method is a wee bit more nuanced, so check the documentation for more detail"""
        signal = self.data_snrpshot
        
        if hasattr(self, 'snr_threhold'):
            threshold = self.snr_threshold
        else:
            self.snr_threhold = self.calculate_snr_threshold()
            threshold = self.snr_threshold
        

        signal_bool = np.full(signal.shape,False,dtype=bool )
        
        for i in range(0,signal.shape[0]):
            for j in range(0,signal.shape[1]):
                
                if signal[i,j] > threshold:
                    signal_bool[i,j] = True
                
                if i > 0 and j > 0:
                    signal_near = signal[i-1:i+2,j-1:j+2]
                    signal_near = signal_near > threshold
                    local_thresh = 6
                    
                elif i == 0 and j > 0:
                    signal_near = signal[i:i+2,j-1:j+2]
                    signal_near = signal_near > threshold
                    local_thresh = 3
                elif i > 0 and j == 0:
                    signal_near = signal[i-1:i+2,j:j+2]
                    signal_near = signal_near > threshold
                    local_thresh = 3
                elif i == 0 and j == 0:
                    signal_near = signal[i:i+2,j:j+2]
                    signal_near = signal_near > threshold
                    local_thresh = 1 
                if signal_near.sum() < local_thresh:
                        signal_bool[i,j] = False
                    
        return signal_bool
    
    def calculate_optimal(self):
        """
        This function calculates and finds the optimal timing pairs for both 4 pulse and 5 pulse DEER.
        """


        if hasattr(self, 'data_snrpshot'):
            norm_signal = self.data_snrpshot
        else:
            try:
                self.snr_normalize()
            except RuntimeError:
                print("Could not auto SNR normlaize the data, please do so manual with the snr_normalize method")
        
        v_snr_bool = self._stability_check()
        time = self.time
        
        # Calculate optimal for 4-pulse
        tau2_index = np.any(v_snr_bool,axis=0).shape[0] - np.flip(np.any(v_snr_bool,axis=0)).argmax() -1
        tau1_index = norm_signal[:,tau2_index].argmax()
        self.index_4p = (tau1_index,tau2_index)
        self.time_4p = (time[0][tau1_index],time[1][tau2_index])

        # Calculate optimal for 5-pulse
        tau1_mesh,tau_2_mesh = np.meshgrid(time[0],time[1])
        tau1ptau2_grid = tau1_mesh + tau_2_mesh
        np.place(tau1ptau2_grid,~v_snr_bool,0)
        max_pos = np.argwhere(tau1ptau2_grid == tau1ptau2_grid.max())
        if max_pos.shape[0] ==1: 
            self.index_5p = (tuple(max_pos[0]))
        else:
            mpT=max_pos.T
            self.index_5p = mpT[:,norm_signal[mpT[0],mpT[1]].argmax()]

        self.time_5p = (time[0][self.index_5p[0]],time[1][self.index_5p[1]])

        return True

    def calculate_snr_threshold(self):
        """Quick script to calculate the SNR threshold from imported values"""
        num_of_acq = self.time_target*3600 /(self.shrt * 1e-6) 
        self.snr_threshold = self.snr_target / np.sqrt(num_of_acq / (self.trace_length))
        return self.snr_threshold 

    def set_snr_threshold(self,value):
        """Sets the SNR_threshold"""
        self.snr_threshold = value

    def set_snr_target(self,value):
        self.snr_target = value

    def set_time_target(self,value):
        self.time_target = value

    def estimated_snr(self,exp):
        num_of_acq = self.time_target*3600 /(self.shrt * 1e-6) 
        if exp == '4-pulse':
            return self.data_snrpshot[self.index_4p] * np.sqrt(num_of_acq / (self.trace_length))
        elif exp == '5-pulse':
            return self.data_snrpshot[self.index_5p] * np.sqrt(num_of_acq / (self.trace_length))
    
    def snr_normalize(self,shots='auto'):
        """
        Calculates the normalized signal. If the total number of shots cannot be calculated through the Bruker metadata
        file then an inputted value is required. If a value is given then this is given preference over any Bruker values.
        This script assumes a 16 step phase cycle was applied for the 2D decoherence experiment.
        Paramaters:
        shots(Int)='auto'': Number of shots per point, if None then total_shots = 1. If shots = None, then it doesn't update
        the main variable

        """

        if shots == None:
            total_shots = 1
        elif shots == 'auto':
            total_shots = self.scans * self.shots * 16
        else:
            total_shots = shots
        self.calculate_noise_level()
        if shots != None:
            self.data_snrpshot = np.real(self.data) / (np.sqrt(total_shots) * self.noise)
        else:
            return np.real(self.data) / self.noise

    def calculate_noise_level(self):
        """This finds the noise level from std of top right hand corner of the 2D experiment. i.e very long delays. If the
        std deviation is not stable enough then a warning will be thrown"""

        #TODO: Add warning and noise level checking
        signal = np.real(self.data)
        width,height = signal.shape
        start_frac = self.noise_frac # This needs to be generated and checked for each data set to make sure it is valid
        self.noise = np.std(signal[int(np.floor(width*start_frac)):,int(np.floor(height*start_frac)):])


    def create_twoD_plot(self,norm=None,contour=None,optimal=False,**kwargs):
        """
        This is the general method for generating 2D plots of these 2D Decoherence plots. There are many optional input parmeters which will dictate how they plot. 
        Some of these possible outputs may require more infomation then is automatically filled from Deerload. If this is a case it will return an exception.

        Parameters:
        norm(string or None): The normalization of the fit. Options: None,'Tau2','SNR','SNRpShot'
        contour(string or None): How the contours are plotted. Options: 'Auto', 'SNRpShot','None'
        optimal(bool): Should the optimal 4pulse and 5pulse pairs be plotted.

        Keywords:
        cmap: To set a particular colormap
        axis: To return a subplot
        title: A plot title
        """


        if not(hasattr(self,'snr_threshold')):
            self.calculate_snr_threshold()

        if contour == 'SNRpShot':
            contour_scale = np.logspace(np.log10(self.snr_threshold),np.log10(self.snr_threshold*64),7,base=10)
            contour_norm = None  
        elif contour == 'Auto':
            contour_scale = 15
            contour_norm = colors.PowerNorm(gamma=0.5)


        
        if norm == 'Tau2':
            signal = np.real(self.data) / np.real(self.data).max(axis=1)[:,None]
            cbar_label = None
            contour_scale = 15
            contour_norm = None
            grid_norm = None
        elif norm == 'SNR':
            signal = self.snr_normalize(shots=None)
            grid_norm = colors.PowerNorm(gamma=0.5)
            cbar_label = 'SNR'
            
        elif norm =='SNRpShot':
            if not(hasattr(self,'data_snrpshot')):
                self.snr_normalize()
            signal = self.data_snrpshot
            cbar_label = r'SNR /$n^{-1/2}$'
            grid_norm = colors.PowerNorm(gamma=0.5)
        else: #No Normalization
            signal = np.real(self.data)
            grid_norm = None
            cbar_label = None
        
        if optimal:
            if not(hasattr(self,'index_4p')):
                self.calculate_optimal()
        
        # Start to generate the figure
        if 'cmap' in kwargs:
            cmap = kwargs['cmap']
        else:
            cmap = cm.get_cmap(name='Reds', lut=None)
        cmap2 = cm.get_cmap(name='gist_gray', lut=None) #Cmap for the contours
        scale = [min(self.time[0]),max(self.time[0]),min(self.time[1]),max(self.time[1])]
        
        if 'figure' in kwargs:
            fig = kwargs['figure']
            axs = fig.subplots(1,1)
        else:
            fig = plt.figure(figsize=(6,6),dpi=150)
            axs = fig.subplots(1,1)
        if contour != None:
            axs.contour(signal,extent=scale,levels=contour_scale,cmap=cmap2,alpha=0.9,zorder = 5,norm=contour_norm)
        im  = axs.pcolormesh(self.time[0],self.time[1],signal,shading='auto',alpha=0.9,cmap=cmap,zorder = 0, norm=grid_norm)
        axs.plot([min(self.time[0]),max(self.time[0])],[min(self.time[0]),max(self.time[0])],color='k',linestyle='--',alpha=0.9,zorder=6)
        axs.set_xlabel(r'$\tau_1 / (us)$')
        axs.set_ylabel(r'$\tau_2 / (us)$')
        axs.set_xlim(min(self.time[0]),max(self.time[0]))
        axs.set_ylim(min(self.time[1]),max(self.time[1]))
        axs.set_aspect('equal')
        if optimal:
            axs.scatter(self.time_4p[0],self.time_4p[1],c='b',s=15,label='4-pulse',zorder=10)
            axs.scatter(self.time_5p[0],self.time_5p[1],c='lime',s=15,label = '5-pulse',zorder=10)
            plt.legend(loc=2)
            axs.text(0,-0.1,fr'$\tau_1 = ${self.time_4p[0]}$\mu s$ & $\tau_2 = ${self.time_4p[1]}$\mu s$',c='b',transform=axs.transAxes)
            axs.text(0.55,-0.1,fr'$\tau_1 = ${self.time_5p[0]}$\mu s$ & $\tau_2 = ${self.time_5p[1]}$\mu s$',c='lime',transform=axs.transAxes)

        
        if 'title' in kwargs:
            axs.set_title(kwargs['title'])
        
    

        cbar = fig.colorbar(im, label =cbar_label )
        
        return fig

    def create_slice_plot(self,axis,target,norm=None,**kwargs):
        """
        This function returns a plot of a slice of the 2D decoherence experiment

        parameters:
        axis(0 or 1): The direction of the the slice (0 = constant tau2; 1 = constant tau1)
        target(float): in us the target time for the slice, due to varying spacing this will find the closest time.
        """


        if norm == 'Tau2':
            signal = np.real(self.data) / np.real(self.data).max(axis=1)[:,None]
            norm_label = None  
        elif norm == 'SNR':
            signal = self.snr_normalize(shots=None)
            norm_label = 'SNR'
    
        elif norm =='SNRpShot':
            if not(hasattr(self,'data_snrpshot')):
                self.snr_normalize()
            signal = self.data_snrpshot
            norm_label = r'SNR /$n^{-1/2}$'
        else: #No Normalization
            signal = np.real(self.data)
            norm_label = None

        fig = plt.figure(figsize=(6,6),dpi=150)
        axs = fig.subplots(1,1)

        slice_pos = abs(self.time[axis]-target).argmin()
        time = self.time[axis][slice_pos]
        if axis == 0: # This means that we are varying tau1 and have a constant tau2
            slice_data = signal[:,slice_pos]
            slice_time = self.time[int(not(bool(axis)))]
            axs.set_xlabel(r'$\tau_1 / (us)$')
        elif axis == 1: # This means that we are varying tau2 and have a constant tau 1
            slice_data = signal[slice_pos,:]
            slice_time = self.time[int(not(bool(axis)))]
            axs.set_xlabel(r'$\tau_2 / (us)$')

        else:
            raise ValueError('Axis must be either 0 or 1')

        axs.plot(slice_time,slice_data)

        if 'title' in kwargs:
            axs.set_title(kwargs['title'])

        return fig , time
    
    
    
    def invert_signal(self):
        """
        This function completely inverts the signal. I.e. -x -> x and -y -> y.

        This is usefull if the measured sameple is pi out of phase, or the phase was aligned for a Hahn echo not a refocused echo
        """

        try:
            self.data = self.data * -1
        except:
            print('Please load data')

            
    def _data_transpose(self):
        """
        This is an internal method that flips the data along the diagnal.
        
        """
        
        self.data  = self.data.T

    def value_at_pos(self,pos:tuple,norm=None):
        """
        Function that finds the value in this 2D experiment at specified coordinates
        Parameters:
        pos(tuple): (tau1 pos, tau2 pos)
        norm=Nome: Can be used to select the normalisation. Options are: None, 'SNRpShot','SNR'  
        """

        if norm == None:
            return self.data(pos[0],pos[1])
        elif norm == 'SNRpShot':
            return self.data_snrpshot(pos[0],pos[1])
        elif norm == 'SNR':
            data = self.snr_normalize(shots=None)
            return data(pos[0],pos[1])
        else:
            raise ValueError('Normalization option value error')

    def value_at_time(self,pos:tuple,norma=None):
        """
        Function that finds the closest value in this 2D experiment for specified delay times
        Parameters:
        pos(tuple): (tau1 time, tau2 time)
        norm=Nome: Can be used to select the normalisation. Options are: None, 'SNRpShot','SNR'  
        """
        # First we need to find the closest posible values
        tau1_pos = np.argmin(abs(self.time[0]-pos[0]))
        tau1_time = self.time[0][tau1_pos]
        tau2_pos = np.argmin(abs(self.time[1]-pos[1]))
        tau2_time = self.time[1][tau2_pos]
        print(f'Closest times: \u03C41 ={tau1_time} & \u03C42 ={tau2_time}')

        return self.value_at_pos((tau1_pos,tau2_pos),norm=norma)
