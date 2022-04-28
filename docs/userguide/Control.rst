Control
====================

The control element of this package exists inside the hardware module. The long term aim is for this to be easily expanable to homebuilt spectrometer, for a variety of reasons 
we are however currently focused on the spectrometers that we have in house (i.e at ETH Zurich).

Since we can not control a Bruker spectrometer directly we utilise the Bruker Xepr API, in effect we are controling Xepr that in turns control the spectrometer. Subsequently, a
Xepr must still be started and connected for this to work. A spectrometer with up to date software is also required to open API gateway. 




Bruker
-------------------

.. warning:: 
    The stardard Bruker Safety Check must be completed before autoDeer is started. AutoDeer can not check for protection switches and glass modes. 

Before a autoDeer can connect to Xepr, the gateway must be first be opened in Xepr. ::

        import autoDeer.hardware.xepr_api_adv as api

        xepr=api()
        xepr.find_Xepr()
        xepr.find_cur_exp()
        xepr.find_hidden()



Hybrid Systems
-------------------

Hybrid systems are often unique so here we present how to connect to a Keysight M8190A AWG and how to send both nutation experiments
and probe pulses to it. This is not designed to complete. ::
    
        awg_ip = '129.132.218.87'

        awg = ksawg.Interface()
        awg.open_con()
        awg._inst_init(awg_ip)
        print(f"Checking Sequencing mode:{awg.getFunctionMode(3)}")

        sampling_freq = 12e9
        sampling_period = 1/sampling_freq
        grad = 64


.. toctree::
    :hidden:
    :caption: Bruker

    ./Bruker