Control
====================

The control element of this package exists inside the hardware module. The long term aim is for this to be easily expanable to homebuilt spectrometer, for a variety of reasons 
we are however currently focused on the spectrometers that we have in house (i.e at ETH Zurich).

Since we can not control a Bruker spectrometer directly we utilise the Bruker Xepr API, in effect we are controling Xepr that in turns control the spectrometer. Subsequently, a
Xepr must still be started and connected for this to work. A spectrometer with up to date software is also required to open API gateway. 