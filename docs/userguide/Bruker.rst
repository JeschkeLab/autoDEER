Bruker
====================


Running an experiment on a Bruker system is the simplest, though least powerful method.
The procedure is built upon PulseSpel and has been designed to be easy to expand to new experiments and customize.

Connecting
-------------------
Connecting to Xepr is a simple process, however it does require that the Bruker API is enabled on Xepr.


::

        import autoDeer.hardware.xepr_api_adv as api

        xepr=api()
        xepr.find_Xepr()
        xepr.find_cur_exp()
        xepr.find_hidden()

Setting key experiment parameters
---------------------------------------------

Whilst, the API can control/Read almost all Xepr features and this is still possible through autoDeer, the most common and useful parameters have had a wrapper
built around them to simplify them these include:

==================       =============   ============
Parameter                 Read             Write
==================       =============   ============ 
Bridge Frequency           X               O         
Frequency Counter              X              O        
Field                    X               X       
==================       =============   ============ 


Running an experiment
-------------------------
Once the field and frequency is set, a pulsed experiment will likely be next. This can be done using the ``api.run_general`` command. This command can be used
to set the PulseSpel files, phase cycling, dimensions, variables and acquisition settings. PulseSpel files for the most common EPR experiments are included,
otherwise users can use their own by specifying the full path. Here is an example. ::

    import autoDeer.hardware.xepr_experiments as exp 

    exp.run_general(
        xepr,                                                       # Bruker API class
        ["/PulseSpel/HUKA_DEER_AWG"],                               # Local Path to PulseSpel files
        ["Field Sweep +<x>/-<x> pg=200ns","Field sweep +<x>/-<x>"], # Name of both experiment and phase cycling
        {"PhaseCycle":False,"ReplaceMode":False},                   # Dictionary containing acquisition settings 
        {"p0":16,"p1":16,"h":25,"n":1,"d0":600}                     # Dictionary of pulseSpel variables to change
        )



**Currently supported acquisition settings**

=============================           ===================   
Acquisition Setting                     Default             
=============================           ===================    
On-Board Phase cycling                  True (on)                        
Replace Mode                            False (off)                      
Acquisition mode                        1 (PulseSpel)                      
=============================           ===================    


Reading a Dataset
-------------------------


Now that an experiment is running, the data will need to be read into Python. This can be done concurrently and at the end. Here we have three different commands:

    1. 'api.acquire_dataset()' - This acquires the data now!
    2. 'api.acquire_scan()' - This acquires the data at the end of the current scan.
    3. 'api.acquire_scan_at(scan_num)' - This acquires the data at a given scan number.
   
