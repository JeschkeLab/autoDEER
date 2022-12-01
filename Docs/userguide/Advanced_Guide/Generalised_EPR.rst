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

