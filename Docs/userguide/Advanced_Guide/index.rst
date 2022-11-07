Advanced User Guide
==========================
An autoDEER experiment is constructed from a collection of modules; 
These can either be Experimental or Data analysis in nature. Advanced users have
the oportunity to develop their own automatic EPR experiments by modifying these
modules and piecing them together in a python script or Jupyter notebook. 


A modular structure
---------------------------


Building a custom experiment
-----------------------------------------

.. caution::
    A basic understanding of python code is essential for building a custom 
    experiment. The greater your understanding the more powerful the code.

The first requirement in any custom experiment, is to connect to the relevent hardware.
Currently, this is limited to just Bruker Xepr but will expand to homebuilt and 
hybrid systems soon.

.. code:: python
    
    from autoDEER.hardware import xepr_api
    Xepr = xepr_api.connect()

When using a non-AWG system, tuning the channels will be the first task. This
has been implemented for Bruker systmes in the :class:`xepr_experiments` class. 

.. code:: python

    tune = exp.MPFUtune(Xepr, "Hahn", 16, 650)
    tune.tune({'+<x>': "R+", '-<x>': "R-"}, tol=1, bounds=[20, 60])

There are stages to this code. The first line builds the class. Here we pass to
the class the api module (here our Xepr), we then tell it what type of experiment
we would like to use for tuning. This can be either a Hahn or refocused echo. 
Next we pass the target :math:`\pi/2` pulse length, and finally the correct d0 
time. The concept of a d0 time is rather unique to Bruker systems and is explained 
in more detail on the dedicated Bruker page. 

The second line is what actually starts the optimisation. Firstly, we pass a 
dictionary of which channels we want to tune and to what phase we want to tune. 
I.e. `'+<x>': "R+"` means to tune the `+<x>` channel to Real echo up. We can 
also pass extra parmaters such as a optimisation tolereance and place additional
bounds. 


Now that the channels are tuned, we can now run an experiment. In this example,
we will choose a simple Echo Detected Field Sweep. This is likely the first
experiment to be measured on a new sample, a detailed list of other experiemnts
is given on the page `Experiments`

.. code:: python

    exp.run_general(
        xepr,
        ["/PulseSpel/HUKA_DEER_AWG"],
        ["Field Sweep +<x>/-<x> pg=200ns", "Field sweep +<x>/-<x>"],
        {"PhaseCycle": False, "ReplaceMode": False},
        {"p0": 16, "p1": 16, "h": 25, "n": 1, "d0": 660}
        )

Once the experiment has run, analysis can then be done. The first step is
getting the data from the spectrometer into autoDEER. This is done through the 
:code:`xepr.acquire_dataset()` command. There are a few other versions of this
command explained in the `control` section. This command will return an instance
of the class :class:`dataset`, this can then be imported into an analysis module.

.. code:: python

    fs_data = xepr.acquire_scan()
    fs = FieldSweep()
    fs.import_from_dataclass(fs_data)
    fs.find_max()

    xepr.set_field(fs.max_field)
    Bc = xepr.get_field()
    fc = xepr.get_counterfreq()

    fs.plot()
    gyro_exp = fs.calc_gyro(xepr.get_counterfreq())
    print(fs.gyro)

Many of the analysis modules have associated methods that perfom some calculation
or generate a specific plot. In this example, we use three methods. First, we 
have the find max value, we then generate a plot and finally we calculate the
gyromagnetic ratio.

.. note::
    The gyromagnetic ratio is given in units of G/GHz, not the more conventionaly
    T/GHz. 

For building more advanced experiments, please look at the avaliable modules as
well as the sections on the specific hardware APIs.


.. toctree::
    :hidden:
    :maxdepth: 2
    :caption: Modules

    ./Experiments.rst
    ./Analysis.rst
    ./Dask.rst

.. toctree::
    :hidden:
    :maxdepth: 2
    :caption: Hardware Control

    ./Bruker.rst

.. toctree::
    :hidden:
    :maxdepth: 2
    :caption: Developer Guide

    ./Logging.rst
