Generalised EPR
===============

Currently, the most common way to describe a pulsed experiment is using 
PulseSpel. Whilst, for many experiments it represents an easy and sufficently
practical method to describe a sequence, it lacks the easy of use that one might
wish for and the ability to easily describe more complex AWG experiemnts.

Here we present, a python based and styled approach to design and create pulsed
EPR experiemnts.

The general structure
----------------------
In this approach, we build up an experimental sequence in terms of pulses (incl. detection events).

Both the sequences and containing pulses are developed from parameters. Each 
parameter is stored alongside its unit and a brief description to reduce 
potential confusion.

Pulses
++++++++++++++++++++
All pulses share a collection of parameters and methods. Special pulses will 
often require additional parameters.

Pulses vs Delays
------------------------
There are two equiavlent but different ways of defining a pulse sequence: delay
focused or pulse focused. Most EPR textbooks are delay focused whilst most hardware
operate (at a fundamental level) in pulse focused. 

In a delay focused sequence, each pulse element has a distinct length (tp) but not a
position. It's position is always relative to when the previous element finishes,
pulse delays are represented by a delay object of some length. This is also how
a Bruker PulseSpel experiment is written.

In a pulse focused sequence, each pulse element has **both** a distinct length (tp) 
and a position (t). Inter-pulse delays are not distinctly represented but are
intfered from the difference in pulse position. Many waveform generators operate
somewhat in this fasion. 

Both approaches are supported here. Automatically any element without a position (t)
is considered to be of delay focused type and vise versa. Converting between approach
is also possible.

Useful Methods
++++++++++++++++++

- Sequence.convert()

- Sequence.isDelayFocused()

- Sequence.isPulseFocused()

- Pulse.isDelayFocused()

- Pulse.isPulseFocused()

Phase Cycles
------------------------
Phase cycles are considered to be a parameter of the Pulse. From this the 
complete sequence based phase cycle can be created.


Adding Progression
-------------------------


PulseSpel Conversion
------------------------
Bruker Elexsys2 based spectrometers mostly require a PulseSpel script and not
a data structure. 

The first step in the conversion is for the sequence to be full delay focused 
and reduced. 
