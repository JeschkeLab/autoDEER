# This file contains the class for running 2D Decoherence experiments
# all code assoicated with processing these experiments are included in
# this file.

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
import deerlab as dl
import os

class TwoD_Experiment:
    """"
    TwoD_Experiment: This is a class for loading and processing 2D Decoherence experiments
    """

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
            self.time,self.data,self.params = dl.deerload(file,full_output=True)
            self.scans = int(self.params.get('DSL').get('recorder').get('NbScansDone'))
            self.shots = int(self.params.get('DSL')['ftEpr']['ShotsPLoop'])
            self.shrt = float(self.params.get('DSL')['ftEpr']['ShotRepTime'].split()[0]) # in us
            self.trace_length = 512 #TODO: This needs to be moved into an input parameter during initialization of this class.
        else:
            self.time,self.data= dl.deerload(file,full_output=False)
            self.params = None
        
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

        v_snr_bool = self._stability_check()
        time = self.time

        if hasattr(self, 'data_snrpshot'):
            norm_signal = self.data_snrpshot
        else:
            try:
                self.snr_normalize()
            except RuntimeError:
                print("Could not auto SNR normlaize the data, please do so manual with the snr_normalize method")

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
        num_of_acq = 2*3600 /(self.shrt * 1e-6) 
        self.snr_threshold = 20 / np.sqrt(num_of_acq / (self.trace_length))
        return self.snr_threshold 

    def set_snr_threshold(self,value):
        """Sets the SNR_threshold"""
        self.snr_threshold = value
    
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
        start_frac = 0.75 # This needs to be generated and checked for each data set to make sure it is valid
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
        """


        if not(hasattr(self,'snr_threshold')):
            self.calculate_snr_threshold()

        if contour == 'SNRpShot':
            contour_scale = np.logspace(np.log10(self.snr_threshold),np.log10(self.snr_threshold*64),7,base=10)
        elif contour == 'Auto':
            contour_scale = 10
        
        if norm == 'Tau2':
          signal = np.real(self.data) / np.real(self.data).max(axis=1)[:,None]
          cbar_label = None  
        elif norm == 'SNR':
            signal = self.snr_normalize(shots=None)
            cbar_label = 'SNR'
            
        elif norm =='SNRpShot':
            if not(hasattr(self,'data_snrpshot')):
                self.snr_normalize()
            signal = self.data_snrpshot
            cbar_label = r'SNR /$n^{-1/2}$'
        else: #No Normalization
            signal = np.real(self.data)
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

        fig = plt.figure(figsize=(6,6),dpi=150)
        axs = fig.subplots(1,1)
        axs.contour(signal,extent=scale,levels=contour_scale,cmap=cmap2,alpha=0.9,zorder = 5)
        im  = axs.pcolormesh(self.time[0],self.time[1],signal,shading='auto',alpha=0.9,cmap=cmap,zorder = 0, norm=colors.PowerNorm(gamma=0.5))
        axs.plot([min(self.time[0]),max(self.time[0])],[min(self.time[0]),max(self.time[0])],color='k',linestyle='--',alpha=0.9,zorder=6)
        axs.set_xlabel(r'$\tau_1 / (us)$')
        axs.set_ylabel(r'$\tau_2 / (us)$')
        axs.set_aspect('equal')
        if optimal:
            axs.scatter(self.time_4p[0],self.time_4p[1],c='b',s=15,label='4-pulse',zorder=10)
            axs.scatter(self.time_5p[0],self.time_5p[1],c='lime',s=15,label = '5-pulse',zorder=10)
            fig.legend(loc=2)
        

        cbar = fig.colorbar(im, label =cbar_label )
        
        return fig
