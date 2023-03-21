DEER Variants
-------------

Over the years many different versions of DEER spectroscopy have been 
developed. Initially, it was 3-pulse and the 4-pulse. With the advent of 
high-speed AWGs the number of sequences exanded even more to include 5-pulse, 
7-pulse and nDEER.

AutoDEER currently supports these variants:

1. 4-pulse
2. 5-pulse
3. 7-pulse
4. nDEER of all of the above

Other variants can be created using the autoEPR package. If you would like them
to be added to autoDEER, please make a pull request on GitHub.

Selecting a version
+++++++++++++++++++

..  code-block:: python
    sequence = DEERSequence(B=B, LO=LO, reptime=reptime, )


References
++++++++++
