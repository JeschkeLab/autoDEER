# autoDeer
A python package designed for the running of automated and optimised DEER/PELDOR experiments for pulsed EPR. 

## Dependencies
1) Numpy
2) Xepr API for python 3
2) Scipy
3) DeerLab (https://jeschkelab.github.io/DeerLab/)

## Troubleshooting
- <Can't find any open instances of Xepr>
    This has two main causes:
    1) Xepr API hasn't be started. Please do "Processing" -> "XeprAPI" -> "Enable Xepr API"
    2) The kernal has recently lost connection, and a background kernal is connected. Please open <htop> in terminal and
    kill the process of the other python kernal

## Acknowledgements 
In developing the Bruker spectrometer control, insipration and help was gained from previous work by:
- Spinach by Ilya Kuprov et. al.
- CustomXepr by San Schott

