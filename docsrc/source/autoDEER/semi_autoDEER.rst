Semi - autoDEER
===============

What is a semi-autoDEER experiment?
-----------------------------------
A semi-autoDEER experiment is an autoDEER experiment run inside a Jupyter 
notebook where each sub-experiment is run individually this allows an experienced
user to intercede and take control when necessary or provide slight modifications
to the script. It has the same form and layout as a fully autoDEER experiment. 

In this example, we show how to run a five-pulse DEER experiment using 
rectangular pulses.

Setup
-----
First, we must import all the necessary packages and the interface. In this
example we are using a custom ETH-AWG interface.

.. code-cell:: python
    :execution-count: 1

    import matplotlib.pyplot as plt
    import numpy as np
    import autodeer as ad
    from autodeer.hardware.ETH_awg import ETH_awg_interface
    from autodeer.sequences import *
    import autodeer.hardware.criteria as criteria
    import time

Here we define a few starting parameters as well as initialise the interface.
It is at this step that a connection between autoDEER and your spectrometer API
is made.

.. code-cell:: python
    :execution-count: 2

    gyro_N = 0.0028087 # This is just an initial guess of the gyromagnetic ratio, a more precise value will be calculated from the field sweep. 
    LO = 34.2172 # The is the central frequency of your resonator
    reptime = 2e3 # This should be 3-5e3 when using 50K
    interface = ETH_awg_interface()

Field Sweep
-----------

An Echo Detected Field Sweep (EDFS) is detected. This allows a more precise 
gyromagnetic ratio to be calculated. Before we can do this, we must first 
tune a pair of rectangular pulses.

.. code-cell:: python
    :execution-count: 3

    p90, p180 = interface.tune_rectpulse(tp=12, LO=LO, B=LO/gyro_N, reptime = reptime)

.. code-cell:: python
    :execution-count: 4

    fsweep = FieldSweepSequence(
        B=LO/gyro_N, LO=LO,reptime=2e3,averages=1,shots=120,
        Bwidth = 500, 
        pi2_pulse=p90, pi_pulse=p180,
        )
    print(fsweep)

.. collapse:: Output

    .. output-cell::
        :execution-count: 4

        ###############################################################################
        autoDEER Sequence Definition
        ###############################################################################
        Sequence Parameters: 
        Name       Value        Unit       Description                    
        B          12183        Gauss      The static B0 field for the experiment 
        LO         34.217       GHz        The local oscillator frequency. 
        reptime    2000         us         The shot repetition time       
        averages   1            None       The number of averages to perform. 
        shots      120          None       The number of shots per scan.  
        det_window 128          None       The length of the default detection gate 
        time       0:4:0        HH:MM:SS   Estimated sequence run time    
        Bwidth     500          Gauss      Field sweep width              

        Events (Pulses, Delays, etc...): 
        iD   t (ns)   tp (ns)  scale    type           Phase Cycle                             
        0        0       12     0.14    RectPulse      [+(+x) -(-x)]                           
        1      500       24     0.14    RectPulse      [+(+x)]                                 
        2     1000      128    N/A      Detection                                              

        Progression: 
        Pulse      Prog. Axis Parameter  Step       Dim        Unit       
        Seq        0          B          1          501        Gauss      
        ###############################################################################
        Built by autoDEER Version: 0.4_dev
        ###############################################################################

.. code-cell:: python
    :execution-count: 5

    interface.launch(fsweep,savename="autoDEER_fieldsweep",IFgain=2)
    while interface.isrunning():
        time.sleep(10)
    dataset = interface.acquire_dataset()
    fsweep_analysis = ad.FieldSweepAnalysis(dataset)
    gyro_exp = fsweep_analysis.calc_gyro(LO)
    fsweep_analysis.plot();

Resonator Profile
-----------------

A resonator profile allows us to find the central frequency and optimise the 
frequency so that most of the spins are inside our 'dip'.  

.. code-cell:: python
    :execution-count: 6

    RPseq = ResonatorProfileSequence(
        B=LO/gyro_exp, LO=LO,reptime=2e3,averages=1,shots=200,
        pi2_pulse=p90, pi_pulse=p180,
    )
    print(RPseq)

.. collapse:: Output

    .. output-cell::
        :execution-count: 6

        ###############################################################################
        autoDEER Sequence Definition
        ###############################################################################
        Sequence Parameters: 
        Name       Value        Unit       Description                    
        B          12182        Gauss      The static B0 field for the experiment 
        LO         34.217       GHz        The local oscillator frequency. 
        reptime    2000         us         The shot repetition time       
        averages   1            None       The number of averages to perform. 
        shots      200          None       The number of shots per scan.  
        det_window 128          None       The length of the default detection gate 
        time       0:13:38      HH:MM:SS   Estimated sequence run time    

        Events (Pulses, Delays, etc...): 
        iD   t (ns)   tp (ns)  scale    type           Phase Cycle                             
        0        0        4        1    RectPulse      [+(+x)]                                 
        1     2000       12     0.14    RectPulse      [+(+x) -(-x)]                           
        2     2500       24     0.14    RectPulse      [+(+x)]                                 
        3     3000      512    N/A      Detection                                              

        Progression: 
        Pulse      Prog. Axis Parameter  Step       Dim        Unit       
        0          0          tp         2          33         ns         
        Seq        1          B          7.1201     31         Gauss      
        Seq        1          LO         0.02       31         GHz        
        ###############################################################################
        Built by autoDEER Version: 0.4_dev
        ###############################################################################

.. code-cell:: python
    :execution-count: 7

    interface.launch(RPseq, savename="autoDEER-resonator_profile", IFgain=2)
    while interface.isrunning():
        time.sleep(10)
    dataset = interface.acquire_dataset()
    respro = ad.ResonatorProfileAnalysis(
        nuts = dataset.data.T,
        freqs = RPseq.LO.prog[0][1],
        dt=2
    )
    respro.process_nutations(threshold=1)
    respro.plot(fieldsweep=fsweep_analysis);

DEER
----

Now we have the necessary setup experiments it is time to start thinking about
DEER. A DEER experiment is normally formed from three types of pulses. An 
excitation π/2 pulse, a refousing π pulse and finally a pump π pulse. There will
be a frequency seperation between the excitation pulse and the pump pulse. 

.. code-cell:: python
    :execution-count: 8

    exc_pulse = RectPulse(  
                    tp=12, freq=0, flipangle=np.pi/2, scale=0
                )
    ref_pulse = RectPulse(  
                    tp=12, freq=0, flipangle=np.pi, scale=0
                )
    pump_pulse = RectPulse(  
                    tp=12, freq=0, flipangle=np.pi, scale=0
                )
    det_event = Detection(tp=512, freq=0)

.. code-cell:: python
    :execution-count: 9

    fpump,fobs = ad.calc_optimal_deer_frqs(
        fsweep_analysis,pump_pulse,exc_pulse,exc_limits=(-0.08,0.025))
    exc_pulse.freq.value += fobs
    ref_pulse.freq.value += fobs
    det_event.freq.value += fobs
    pump_pulse.freq.value += fpump
    ad.plot_optimal_deer_frqs(fsweep_analysis,pump_pulse,exc_pulse);


.. code-cell:: python
    :execution-count: 10

    p90, p180 = interface.tune_rectpulse(
        tp =12, LO=LO+exc_pulse.freq.value, B=LO/gyro_exp, reptime = reptime)
    p90.freq.value = exc_pulse.freq.value
    p180.freq.value = p90.freq.value
    CPseq = CarrPurcellSequence(
        B=LO/gyro_exp, LO=LO,reptime=2e3,averages=1,shots=200,
        n=2,tau=32e3, pi2_pulse=p90, pi_pulse=p180, det_event=det_event)
    print(CPseq)
    interface.launch(CPseq, savename="autoDEER-CarrPurcell", IFgain=2)

.. collapse:: Output

    .. output-cell::
        :execution-count: 10

        ###############################################################################
        autoDEER Sequence Definition
        ###############################################################################
        Sequence Parameters: 
        Name       Value        Unit       Description                    
        B          12182        Gauss      The static B0 field for the experiment 
        LO         34.217       GHz        The local oscillator frequency. 
        reptime    2000         us         The shot repetition time       
        averages   1            None       The number of averages to perform. 
        shots      200          None       The number of shots per scan.  
        det_window 128          None       The length of the default detection gate 
        tau        32000        ns         Total sequence length          
        n          2            None       The number of pi pulses        
        time       0:10:16      HH:MM:SS   Estimated sequence run time    

        Events (Pulses, Delays, etc...): 
        iD   t (ns)   tp (ns)  scale    type           Phase Cycle                             
        0        0       12     0.16    RectPulse      [+(+x) -(-x)]                           
        1     8000       24     0.16    RectPulse      [+(+x)]                                 
        2    24000       24     0.16    RectPulse      [+(+x)]                                 
        3    32000      512    N/A      Detection                                              

        Progression: 
        Pulse      Prog. Axis Parameter  Step       Dim        Unit       
        1          0          t          10         770        ns         
        2          0          t          30         770        ns         
        3          0          t          40         770        ns         
        ###############################################################################
        Built by autoDEER Version: 0.4_dev
        ###############################################################################

.. code-cell:: python
    :execution-count: 11

    while interface.isrunning():
        time.sleep(10)
    dataset = interface.acquire_dataset()
    dataset.axes[0] = dataset.axes[0]/1e3
    CP_data = ad.CarrPurcellAnalysis(dataset)
    CP_data.fit()
    CP_data.plot();
    tau = CP_data.find_optimal(
        averages=800,SNR_target=20/0.3, target_time=4,
        target_shrt=reptime*1e-6, target_step=0.015)
    print(tau)

.. code-cell:: python
    :execution-count: 12

    exc_pulse = p90
    ref_pulse = interface.tune_pulse(
        ref_pulse, mode="amp_nut", B=LO/gyro_exp, LO=LO, reptime=reptime )
    pump_pulse = interface.tune_pulse(
        pump_pulse, mode="amp_nut", B=LO/gyro_exp, LO=LO, reptime=reptime )


.. code-cell:: python
    :execution-count: 13

    deer = DEERSequence(
        B=LO/gyro_exp, LO=LO,reptime=2e3,averages=30,shots=80,
        tau1=tau, tau2=tau, tau3=0.2, dt=15,
        exc_pulse=exc_pulse, ref_pulse=ref_pulse,
        pump_pulse=pump_pulse, det_event=det_event)
    deer.five_pulse()
    deer.select_pcyc("16step_5p")
    print(deer)

.. collapse:: Output

    .. output-cell::
        :execution-count: 13

        ###############################################################################
        autoDEER Sequence Definition
        ###############################################################################
        Sequence Parameters: 
        Name       Value        Unit       Description                    
        B          12182        Gauss      The static B0 field for the experiment 
        LO         34.217       GHz        The local oscillator frequency. 
        reptime    2000         us         The shot repetition time       
        averages   30           None       The number of averages to perform. 
        shots      80           None       The number of shots per scan.  
        det_window 128          None       The length of the default detection gate 
        time       6:49:36      HH:MM:SS   Estimated sequence run time    

        Events (Pulses, Delays, etc...): 
        iD   t (ns)   tp (ns)  scale    type           Phase Cycle                             
        0        0       12     0.16    RectPulse      [+(+x)]                                 
        1     2400       12      0.3    RectPulse      [+(+x)]                                 
        2     2600       12     0.36    RectPulse      [+(+x) -(+y) +(-x) -(-y)]               
        3     2800       12      0.3    RectPulse      [+(+x) +(+y) +(-x) +(-y)]               
        4     7800       12     0.36    RectPulse      [+(+x)]                                 
        5    10400      512    N/A      Detection                                              

        Progression: 
        Pulse      Prog. Axis Parameter  Step       Dim        Unit       
        3          0          t          15         320        ns         
        ###############################################################################
        Built by autoDEER Version: 0.4_dev
        ###############################################################################


.. code-cell:: python
    :execution-count: 14

    interface.launch(deer,savename="autoDEER_4hr",IFgain=2)
    interface.terminate_at(
        criteria.DEERCriteria(tau1=tau,tau2=tau,tau3=0.2,mode="speed"),verbosity=2)

.. code-cell:: python
    :execution-count: 15

    while interface.isrunning():
        time.sleep(10)

    dataset = interface.acquire_dataset()
    fit,ROI,new_tau = ad.DEERanalysis(
        dataset.axes[0]/1000 - tau, dataset.
        data, tau, tau, tau3=0.2,
        num_points=100, compactness=True, plot=False, precision="speed")
    ad.DEERanalysis_plot(fit,background=False,ROI=ROI);

.. code-cell:: python
    :execution-count: 16

    tau = np.around(new_tau/2,1)
    deer = DEERSequence(
        B=LO/gyro_exp, LO=LO,reptime=2e3,averages=1000,shots=50,
        tau1=tau*1e3, tau2=tau*1e3, tau3=0.2*1e3, dt=15,
        exc_pulse=exc_pulse, ref_pulse=ref_pulse,
        pump_pulse=pump_pulse, det_event=det_event)
    deer.five_pulse()
    deer.select_pcyc("16step_5p")
    print(deer)

.. collapse:: Output

    .. output-cell::
        :execution-count: 16

        ###############################################################################
        autoDEER Sequence Definition
        ###############################################################################
        Sequence Parameters: 
        Name       Value        Unit       Description                    
        B          12181        Gauss      The static B0 field for the experiment 
        LO         34.217       GHz        The local oscillator frequency. 
        reptime    2000         us         The shot repetition time       
        averages   1000         None       The number of averages to perform. 
        shots      50           None       The number of shots per scan.  
        det_window 128          None       The length of the default detection gate 
        time       332:0:0      HH:MM:SS   Estimated sequence run time    

        Events (Pulses, Delays, etc...): 
        iD   t (ns)   tp (ns)  scale    type           Phase Cycle                             
        0        0       12     0.14    RectPulse      [+(+x)]                                 
        1     5600       12      0.3    RectPulse      [+(+x)]                                 
        2     5800       12     0.32    RectPulse      [+(+x) -(+y) +(-x) -(-y)]               
        3     6000       12      0.3    RectPulse      [+(+x) +(+y) +(-x) +(-y)]               
        4    17400       12     0.32    RectPulse      [+(+x)]                                 
        5    23200      512    N/A      Detection                                              

        Progression: 
        Pulse      Prog. Axis Parameter  Step       Dim        Unit       
        3          0          t          15         747        ns         
        ###############################################################################
        Built by autoDEER Version: 0.4_dev
        ###############################################################################