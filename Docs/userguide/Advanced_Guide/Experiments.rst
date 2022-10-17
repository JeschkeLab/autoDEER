Experiments
================

autoDEER is capable of running many different EPR experiments. Most of these 
are necessary for conducting and setting up an automatic DEER measurement.
All of these experiments are accessible to the user in Jupyter, meaning that 
other experiments can be automated and run using this package. 


Field Sweep
--------------

In nearly all cases the first measurement to be recorded on a pulse spectrometer
is an Echo Detected Field Sweep (EDFS). This is the zeroth harmonic EPR spectrum and is
equivalent to the integral of Continuous Wave (CW) spectrum. 

The field sweep can also be shown on an approximate frequency axis, allowing
for it to be directly compared with the resonator profile.

A field sweep is formed from a simple two pulse Hahn echo at a constant frequency,
the magnetic B_0 field is then swept.


Resonator Profile
-------------------
Another equally important but, less measured experiment is the resonator profile.
Just like it says on the tin, this is the frequency profile of the resonator and can
be used to accurately determine the central frequency and optimise the field such
that the spectum is located inside the profile("dip").


Pulse Profiles
^^^^^^^^^^^^^^^^^
Pulse profiles are very similar to a resonator profile except instead of giving
infomation about the resonator they can be used to show the behaviour of a test
pulse in the frequency domain. 

Unlike a resonator profile the test pulse isn't nutated or frequency swept.
Just the two pulses that form the Hahn echo have their frequency steped. 

When using a BRUKER (non-AWG) spectrometer this is only possible for the ELDOR 
channel, as the frequency can not be adjusted for individual pulses.


Nutation Experiment
----------------------
At the heart of a Resonator Profile is the nutation experiment however, this
is often a highly useful experiment in its own right. Especially,
when considering incoherent pulses. 

A nutation experiment is formed from a test pulse followed by a Hahn echo. The
test pulse is increased in length whilst all other pulses remain constant.
When the test pulse is of zero length, there will be a maximal signal amplitude
as this is just the standard Hahn echo situation, however as it increases in
length the spin magnitisation has increasingly rotated before the Hahn echo
starts. For example, when the test pulse is at the correct pi/2 length the 
initial magnetisation for the Hahn echo is alrady fully in My, and the whole 
sequence is effectively just 2 pi pulses, so no net transverse magnetisation
can be detected. Equally, when the test pulse is at the correct pi length, the
initial magnetisation for the Hahn ehco is -Mz, so a echo is refocused with
negative amplitude, as it is refocused along -My. 




Two-Pulse ESEEM
------------------



Two Dimensional
^^^^^^^^^^^^^^^^^^^



Auto Tune
------------------


MPFU
^^^^^^^^^^^^


ELDOR
^^^^^^^^^^^^
On a BRUKER spectrometer the ELDOR channel is incoherent with respect to the 
digitiser. This means it is not possible to refocus a hahn echo
with it. Subsequently, the only method to tune the signal is through a 
nutation experiment. 



Custom Experiments
--------------------------
Each of these standard experiemnts are only really and simple and easy to use
alias for a custom experiment behind the scenes. As an advanced user you might 
want to implement your own experiments or your own versions of these standard 
experiments. This software package has been designed with this in mind. 

How a custom experiment is run will be very different depending on the hardware
interface used. Here I present a small example using the Bruker Xepr interface.
When using a Bruker interface, these stil exists PulseSpel files at its core,
it is these files that set out the experiment. Whilst these PulseSpel files
cannot be easily modified with code the key parameters can be. For running
automated experiments this is most important.

.. code:: python

    exp.run_general(
        xepr,
        ["/PulseSpel/HUKA_DEER_AWG"],
        ["Field Sweep +<x>/-<x> pg=200ns", "Field sweep +<x>/-<x>"],
        {"PhaseCycle": False, "ReplaceMode": False},
        {"p0": 16, "p1": 16, "h": 25, "n": 1, "d0": 660}
        )