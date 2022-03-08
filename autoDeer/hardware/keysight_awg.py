"""
This is selection of functions and comands for the control of the keysight M8120 AWG
"""
import pyvisa

class interface:

    def __init__(self,protocol) -> None:
        
        if protocol == "VXI":
            self.ResString = "TCPIP0::localhost::inst0::INSTR"
        elif protocol == "HiSLIP":
            self.ResString = "TCPIP0::localhost::hislip0::INSTR"
        elif protocol == "Socket":
            self.ResString = "TCPIP0::localhost::5025::SOCKET"
        else:
            raise ValueError("Incorect Protocol type. Three LAN Protocols are supported:\n"+
            "\"VXI\" : VXI-11 \n"+
            "\"HiSLIP\" : HiSLIP (recomended) This program requires Agilent I/O Libaries Suite to be installed\n"+
            "\"Socket\" : This must be used if I/O libary doesn't support HiSLIP")
        pass

    def open_con(self,name:str) -> None:
        rm = pyvisa.ResourceManager()
        my_instrument = rm.open_resource(name)
        pass


