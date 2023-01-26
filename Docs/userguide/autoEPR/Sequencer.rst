Sequencer
==============

autoEPR can be used to script and develop EPR pulse sequences in Python. 

Pulse sequences are built up of two main object: :ref:`Sequences <Sequence>`
and :ref:`Pulses <Pulse>`. In its simplist form a pulse sequence is just a
train of pulses. Here, we consider a pulse in its general form. I.e. Pulses are
not just MW pulses but also delays and detection events. 

A sequence also contains additional infomation about the experiment, such as:
the LO frequency; the B0 field; Averages and shots; and the shot repetition 
time. 

There are two ways of thinking about and writing down a pulse sequence:

1. Delay focused (Like PulseSpel)
2. Pulse focused (Like Tables)

When many EPR spectrometrists think about pulse sequences they think of pulses
seperated by delays, and these delays increase and decrease in length to change
the position of the pulses. However, most AWGs and spectrometers think in terms
of pulses with a time position and it is this time positon that moves. Both of
these approaches are physically equivalent, and such both approaches can be 
used in autoEPR. However, most interfaces can only take sequences of one type.
Sequences can be converted between each other with the `Sequence.convert()` 
command. 

Creating your first sequence
--------------------------------------
Here I will demonstrate how to build a simple Hahn Echo sequence using 
rectangular pulses.

>>> import autodeer as ad
>>> import numpy as np
>>> from autodeer.openepr import Sequence, RectPulse, Detection
>>> HahnEcho = Sequence(
>>>     name="Hahn Echo", B=12220, LO=34.0,reptime=3e3,
>>>     averages=1, shots=100
>>>     )
>>> HahnEcho.addPulse(RectPulse(t=0, tp=16, freq=0, flipangle=np.pi/2))
>>> HahnEcho.addPulse(RectPulse(t=200, tp=32, freq=0, flipangle=np.pi))
>>> HahnEcho.addPulse(Detection(t=400, tp=32))

This sequece can also be created using delays.

>>> import autodeer as ad
>>> import numpy as np
>>> from autodeer.openepr import Sequence, RectPulse, Delay, Detection
>>> HahnEcho = Sequence(
>>>     name="Hahn Echo", B=12220, LO=34.0,reptime=3e3,
>>>     averages=1, shots=100
>>>     )
>>> HahnEcho.addPulse(RectPulse(tp=16, freq=0, flipangle=np.pi/2))
>>> HahnEcho.addPulse(Delay(tp=200))
>>> HahnEcho.addPulse(RectPulse(tp=32, freq=0, flipangle=np.pi))
>>> HahnEcho.addPulse(Delay(tp=200))
>>> HahnEcho.addPulse(Detection(tp=32))

**Adding Progression**

A single static echo is only of limited use. To make a pulse move we need to 
add a progressive axis. In this example we will increase the delays of the
above sequence from 200ns to 6us in 100ns steps.

>>> HahnEcho.addPulsesProg(
>>>     pulses=[1,3],
>>>     variables=["tp","tp"]
>>>     axis_id = 0,
>>>     axis = np.arange(2e2,6e3,100)
>>>     multipliers = [1,1]
>>>     )

**Adding Phase Cycles**

In nearly all EPR pulse sequences a phase cycle is required to remove 
additional unwanted echoes and coherence pathways. This is possible in
autoEPR without having to multiply out the full cycle, by adding the phase
cycle indivdually to each pulse. 

This is done by adding an extra dictionary to the pulse delceration. This 
dictionary consitsts of two list: phases ("phases"), detection signs("dets").

For example if we wanted to add the standard DC offset cycle on the first pi/2 
pulse we can simply do.

>>> HahnEcho.addPulse(RectPulse(
>>>     tp=16, freq=0, flipangle=np.pi/2,
>>>     pcyc={"phases":[0, np.pi], "dets":[1, -1]}))

**Viewing the sequence**   

The simplist way to look at a pulse sequence is to print it. Printing a 
Sequence will generate a table containing all the import paramaters and 
infomation.

>>> print(HahnEcho)

**Standard Sequences**   

Since a Hahn Echo is such a standard sequence it comes pre-built as
`HahnEchoSequence`. A full list of standard sequences can be found 
:ref:`here <standard_sequences>`


Using Shaped and Chirped pulses
--------------------------------------
So far we have only looked at monochromatic rectangular pulses. It is often
advantages to use wideband shaped and chirped pulses. Chirped pulses allows for
a wider excitation bandwidth at a lower power level, and shaped pulse are used
to reduce excitation side-lobes. 

Many standard pulses are included, however custom pulses can also be created.
A full list of standard pulses can be found 
:ref:`here <standard_sequences>`


Chirped pulses
++++++++++++++++++++
The simplist type of chirped pulse is a linear chip and it is implemented in
the  :ref:`ChirpPulse <ChirpPulse>` method. ChirpPulse requires the paramaters
as RectPulse when being declared as well as two frequency paramaters. Options:
Bandwidth: `BW`; Initial Frequency `init_freq`; Final Frequency `final_freq`.


>>> from autodeer.openepr import *
>>> test_pulse = ChirpPulse(
>>>     tp = 128, init_freq = -0.1, BW=0.2,
>>>     flipangle = np.pi
>>>     )

Just like with Sequences, Pulses can also be printed to get a table of the 
import infomation.

>>> print(test_pulse)

Custom pulses
++++++++++++++++++++
Unlike when making custom sequences it is not recommended that you use an
instance of the class, instead it is recommended that create a new class which
inherits from the class :ref:`Pulse <Pulse>`.

All pulses, have the attributes `AM` and `FM` which describe their shape and 
frequency dependecies. The recommended way of creating these is to create a 
function called `func` and builds `AM` and `FM` for this specific pulse length.
For a rectangular pulse this is simply:

..  code-block:: python

    def func(self, ax):
        nx = ax.shape[0]
        AM = np.ones(nx)
        FM = np.zeros(nx)
        return AM, FM


..  code-block:: python

    from autodeer.openepr import pulse
    class new_pulse(pulse):

        # First we need to create the class init statment
        def __init__(
            self, tp, freq, t=None, flipangle=None, pcyc=[0],
            name=None) -> None:
        # Next we need to fill the parent-class with its infomation.
        Pulse.__init__(
            self, tp=tp, t=t, flipangle=flipangle, pcyc=pcyc, name=name)
        # All non standard parameters are now declared
        self.freq = Parameter("freq", freq, "GHz", "Frequency of the Pulse")
        self.Progression = False
        # Here we build the FM and AM using our custom function.
        self._buildFMAM(self.func)
        pass
        
        # Here we define the custom function.
        def func(self, ax):
            nx = ax.shape[0]
            AM = np.ones(nx)
            FM = np.zeros(nx)
            return AM, FM





What is not currently supported?
--------------------------------------
1. RF pulses (No ENDOR or DNP)

