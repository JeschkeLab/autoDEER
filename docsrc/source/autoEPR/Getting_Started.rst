Getting Started
===================

autoEPR allows the user to script together multiple EPR sequences and
measurements along with active-data processing to develop completely automatic
push-button style experiments.

A simple automated EPR experiment would consist of these main elements:

1. Defining the Sequence
2. Pulse Tuning
3. Starting the experimental sequence
4. Stopping the experimental sequence
5. Data processing and analysis

A more advanced experiment might contain more than one sequence and/or conduct
active data analysis whilst the sequence is still running. 

Setting up your script
--------------------------
autoEPR is a python package. Before use, it must be imported along with the 
correct interface.

..  code-block:: python

    import autodeer
    from autodeer.hardware import BrukerAWG as interface

This imports the whole autoEPR which can be accessed with 
`autoEPR.[class/function]`. Meanwhile we have loaded the 
:ref:`BrukerAWG <BrukerAWG>` interface which can simply be accesed with
`interface`. This is done to keep scripts as general between interfaces.


1. Defining the sequence
++++++++++++++++++++++++++++
The first step is creating a pulse sequence to be run. autoEPR has a built-in
python-based pulse-sequencer. This can be used to easily develop a general 
spectrometer-independent pulse sequence. There is no need to interface with 
Domain Specific Languages (DSL) such as PulseSpel.

Whilst you can always write your own, many common sequences come as default in
autoEPR. For a full list of standard sequences sees the 
:doc:`API Reference <Reference>`. Here we will look at a simple Hahn-Echo 
experiment. 

.. code-block:: python

    sequence = autodeer.HahnEchoSequence(
        B = 12220, LO=34.01, reptime=3e3, averages=1e3, shots=100
    )

    # details about the sequence can be seen by pritning it
    print(sequence)

2. Tuning for the sequence
++++++++++++++++++++++++++++
Before the sequence can be run additional tuning information must be added. This
comes in the form of adding a scale value to each pulse in the sequence. There
are a few ways of tuning a pulse in EPR. autoEPR supports a selection of these
, however, not all spectrometer support all approaches and not all methods work 
for every scenario. 

1. Echo-amplitude (Amplitude sweep)
2. Echo-amplitude (Pulse duration sweep) [Coming soon]
3. Nutation (Amplitude sweep) 
4. Nutation (Pulse duration sweep) [Coming soon]
5. Pulse profiles. [Coming soon]

In this basic example, we will tune using the first method. On Bruker 
spectrometers amplitude sweeps are not available natively however, autoEPR
solves this by conducting multiple measurements at different amplitudes.

.. code-block:: python

    pulse = interface.tune_pulse(sequence.pulses[0], B, LO, reptime, mode="amp_nut")

This returns a modified pulse that now contains the correct scale value.

If only simple rectangular pulses of equal power are required for a setup 
experiment this can be done using the `tune_rect` function. These pulses are
then cached in `interface.pulses` so they can be easily accessed at a later
time. A set of pi and pi/2 pulses in the cache are required for `tune_pulse`,
otherwise, they will be generated.

.. code-block:: python

    pi2_pulse, pi_pulse = interface.tune_rect(tp, B, LO, reptime)

The cache (`interface.pulses`) is a dictionary. Rectangular pulses added to 
the dictionary will have keys in the form of "p[90/180]_[tp]_[freq]" e.g 
"p90_12_34.1". Remember that if the sample is changed or the resonator
adjusted the cache needs to be cleared: `interface.clear_cache()`.

3. Starting the sequence
++++++++++++++++++++++++++++
Now that spectrometer is tuned and the sequence has the appropriate scaling
information we can now launch the actual sequence. Only one piece of information
is needed, a file name.
.. code-block:: python
    
    interface.launch(sequence, savename="My First Sequence")

Acquiring data from the spectrometer is trivally simple. 

.. code-block:: python

    dataset = interface.acquire_dataset()

This returns a :ref:`dataset<dataset>` object that contains not only the `data`,
`axes` and `sequence` but also any extra paramaters from the spectrometer.


4. Stopping the sequence
++++++++++++++++++++++++++++
Deciding when to stop a measurement is much more involved than when to start. 
Normally, this is done by either time or a set number of scans. Both of these
are normally decided before the measurement has started and before the quality
of the sample is truly known. autoEPR allows for more accurate criteria to be
set such as Signal to Noise(SNR) or for dynamic criteria to be used. An 
example of such dynamical criteria would be in autoDEER where the distance 
distribution after data processing is used.

Setting termination criteria is trivially simple with autoEPR. Here is an 
example where we are using SNR criteria of 30.

.. code-block:: python

    criteria = autodeer.SNRCriteria(threshold=30)
    interface.terminate_at(criteria)
    
This will test the measurement every 10 minutes to see if the threshold has been
reached and then it will terminate the experiment. If a more regular or less 
regular interval is required, this can be specified by passing the argument 
`test_interval=` with the time period given in minutes.

If you want to terminate the experiment immediately, this is also possible.
.. code-block:: python

    interface.terminate()


5. Data analysis
++++++++++++++++++++++++++++

Now we have data we need to analyse it. For standard experiments this is simple
as autoEPR contains built in analysis modules however, you are free to process
the data as required.

.. code-block:: python

    dataset = interface.acquire_dataset()
    analysis = CarrPurcellAnalysis(dataset)
    analysis.fit()
    analysis.plot()

Next Steps
--------------------------

More complicated multi-step experiments can be created by combining mutliple 
pulse sequence into one script.

For more infomation check here:

- :doc:`Sequences and Pulses<Sequencer>`
- :doc:`Interfaces<Interface>`
- :doc:`Standard Sequences, Pulses, Interfaces and Analysis<Reference>`
 