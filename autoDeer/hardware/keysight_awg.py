"""
This is selection of functions and comands for the control of the keysight M8120 AWG through a TCP Socket
"""
import re
import socket
from unittest import result

import logging

hw_log = logging.getLogger('hardware.AWG')

class interface:

    def __init__(self,protocol) -> None:
        pass

    def open_con(self,host,port=5025) -> None:
        try:
            self.instr_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        except socket.error as mg:
            print('Failed to create socket. Error code: ' + str(msg[0]) + ' , Error message : ' + msg[1])
            self.instr_socket.close()
        
        print("Socket Created")

        self.wait_time = 0.2 # A delay time variable in seconds
        
        def is_ip(host):
            try:
                socket.inet_aton(host)
            except socket.error:
                return False 
            else:
                return True
        
        if not is_ip(host):
            try:
                self.ip = socket.gethostbyname(host)
            except socket.gaierror:
                print('Hostname could not be resolved')
                self.instr_socket.close()
        else:
            self.ip = host
        
        self.port = port
        self.instr_socket.connect((self.ip,self.port))

    def _inst_init(self):

        try:
            self.instr_socket.sendall(b"*CLS\n")

            self.instr_socket.sendall(b"*RST;*OPC?\n")
            result = self.instr_socket.recv(8)

            self.instr_socket.sendall(b"*IDN?\n")
            idn_result = self.instr_socket.recv(255)
            print(f"Identification String: {idn_result}")

        except socket.error:
            print("Instrument initalisation failed")
            self.instr_socket.close()
        pass

    def _sendSCPI(self,cmd) -> None:

        if type(cmd) is str:
            cmd_b = cmd.encode()
            cmd_s = cmd
        elif type(cmd) is bytes:
            cmd_b = cmd
            cmd_s = cmd.decode()
        
        self.instr_socket.sendall(cmd_b)
        hw_log.info(f'SCPI_Cmd:{cmd_s}')
        
        pass

    def _recvSCPI(self,buffer):
        data_b = self.instr_socket.recv(buffer)
        data = data_b.decode
        hw_log.info(f'SCPI_recv:{data}')
        
        return  data

    def _sendrecvSCPI(self,cmd:str,buffer:int):

        self._sendSCPI(self,cmd)
        data = self._recvSCPI(self,buffer)
        
        return  data


    def Abort(self,chn:int = 0)->None:
        """Abort Stop signal generation on either one or more channels. If channels are coupled bith channels are stopped.

        Parameters
        ----------
        chn : int
            The channel number. 1 or 2 for the respective channel, 3 for both and 0 for default/coupled, by default =0
        """        
        if chn == 3:
            self._sendSCPI(b":ABOR1\n")
            self._sendSCPI(b":ABOR2\n")
        elif chn == 0:
            self._sendSCPI(b":ABOR\n")
        elif chn <3:
            self.__sendSCPI(b":ABOR"+ str(chn).encode + "\n")

    def StartSigGen(self,chn:int = 0)->None:
        """StartSigGen Start signal generation on a channel. Of channels are coupled, both channels are started

        Parameters
        ----------
        chn : int
            The channel number. 1 or 2 for the respective channel, 3 for both and 0 for default/coupled, by default =0
        """
        if chn == 3:
            self._sendSCPI(b":INIT:IMM1\n")
            self._sendSCPI(b":INIT:IMM2\n")
        elif chn == 0:
            self._sendSCPI(b":INIT:IMM\n")
        elif chn <3:
            self.__sendSCPI(b":INIT:IMM1"+ str(chn).encode + "\n")
        
        pass
    
    def getTrigDelay(self,chn:int = 0) -> float:
        """getTrigDelay Get the trigger to data out delay.

        Parameters
        ----------
        chn : int, optional
            The channel number. 1 or 2 for the respective channel and 0 for default/coupled, by default =0

        Returns
        -------
        float
            Read the trigger delay in seconds
        """        
        cmd = b":ARM:DEL"
        buffer = 32
        if chn == 0:
            response = self._sendrecvSCPI(cmd + b"?" +  b"\n",buffer)
        elif chn <3:
            response = self._sendrecvSCPI(cmd + str(chn).encode + b"?" + "\n",buffer)

        return float(response.decode())

    def setTrigDelay(self,delay:float,chn:int = 0) -> float:
        """setTrigDelay Set the trigger to data out delay.

        Parameters
        ----------
        delay : float
            The trigger delay in seconds
        chn : int, optional
            The channel number. 1 or 2 for the respective channel, 3 for both and 0 for default/coupled, by default =0

        Returns
        -------
        float
            Read the trigger delay in seconds
        """

        cmd = b":ARM:DEL"
        buffer = 32
        set_str = str(delay).encode()

        if chn == 3:
            self._sendSCPI(cmd + b"1 " + set_str +  b"\n",buffer)
            self._sendSCPI(cmd + b"2 " + set_str +  b"\n",buffer)

        if chn == 0:
            self._sendSCPI(cmd + b" " + set_str +  b"\n",buffer)
        elif chn <3:
            self._sendSCPI(cmd + str(chn).encode + b" " + set_str+ "\n",buffer)
        

        return self.getTrigDelay(chn)

    def getArmMode(self,chn:int = 0) -> str:
        """getArmMode Get the arming mode

        Parameters
        ----------
        chn : int, optional
            The channel number. 1 or 2 for the respective channel and 0 for default/coupled, by default =0

        Returns
        -------
        str
            Read the arming mode
        """        
        cmd = lambda c: b":INIT:CONT" + str(c).encode() + b":ENABI"
        buffer = 32
        if chn == 0:
            response = self._sendrecvSCPI(cmd("") + b"?" +  b"\n",buffer)
        elif chn <3:
            response = self._sendrecvSCPI(cmd(chn) + str(chn).encode + b"?" + "\n",buffer)

        return response.decode()

    def setArmMode(self,mode:str,chn:int = 0) -> str:
        """setArmMode Set the arming mode

        Parameters
        ----------
        mode : str
            The arming mode. "SELF" or "ARM"ed
        chn : int, optional
            The channel number. 1 or 2 for the respective channel, 3 for both and 0 for default/coupled, by default =0

        Returns
        -------
        str
            Returns the state of the arming mode
        """        
        
        cmd = lambda c: b":INIT:CONT" + str(c).encode() + b":ENABI"
        buffer = 32
        set_str = mode.encode()

        if chn == 0:
            self._sendSCPI(cmd("") + b" " + set_str +  b"\n")
        elif chn ==3:
            self._sendSCPI(cmd(1) + b" " + set_str +  b"\n")
            self._sendSCPI(cmd(2) + b" " + set_str +  b"\n")
            
        elif chn <3:
            self._sendSCPI(cmd(chn) + str(chn).encode + b" " + set_str+ "\n")

        return self.getArmMode(chn)

    def getContMode(self,chn:int = 0) -> int:
        """getContMode Get the continous mode.

        Parameters
        ----------
        chn : int, optional
            The channel number. 1 or 2 for the respective channel and 0 for default/coupled, by default =0

        Returns
        -------
        int
            The continous mode state.
        """        

        cmd = lambda c: b":INIT:CONT" + str(c).encode()
        buffer = 32
        if chn == 0:
            response = self._sendrecvSCPI(cmd("") + b"?" +  b"\n",buffer)
        elif chn <3:
            response = self._sendrecvSCPI(cmd(chn) + str(chn).encode + b"?" + "\n",buffer)

        return int(response.decode())


    def setContMode(self,state:int,chn:int = 0) -> int:
        """getContMode Set the continous mode. This must be used in conjuction with "setGateMode".

        Parameters
        ----------
        state : int
            The state of the continous mode. 0 for off, 1 for on.
        chn : int, optional
            The channel number. 1 or 2 for the respective channel, 3 for both and 0 for default/coupled, by default =0

        Returns
        -------
        int
            The continous mode state.
        """        
        cmd = lambda c: b":INIT:CONT" + str(c).encode() 
        buffer = 32
        set_str = state.encode()

        if chn == 0:
            self._sendSCPI(cmd("") + b" " + set_str +  b"\n")
        elif chn == 3:
            self._sendSCPI(cmd(1) + b" " + set_str +  b"\n")
            self._sendSCPI(cmd(2) + b" " + set_str +  b"\n")

        elif chn <3:
            self._sendSCPI(cmd(chn) + str(chn).encode + b" " + set_str+ "\n")

        return self.getContMode(chn)

    def getGateMode(self,chn:int = 0) -> int:
        """getGateMode Get the continous mode.

        Parameters
        ----------
        chn : int, optional
            The channel number. 1 or 2 for the respective channel and 0 for default/coupled, by default =0

        Returns
        -------
        int
            The Gate mode.
        """        

        cmd = lambda c: b":INIT:GATE" + str(c).encode()
        buffer = 32
        if chn == 0:
            response = self._sendrecvSCPI(cmd("") + b"?" +  b"\n",buffer)
        elif chn <3:
            response = self._sendrecvSCPI(cmd(chn) + str(chn).encode + b"?" + "\n",buffer)

        return int(response.decode())


    def setGateMode(self,state:int,chn:int = 0) -> int:
        """getContMode Set the continous mode. This must be used in conjuction with "setContMode".

        Parameters
        ----------
        state : int
            The state of the gate mode. 0 for off, 1 for on.
        chn : int, optional
            The channel number. 1 or 2 for the respective channel, 3 for both and 0 for default/coupled, by default =0

        Returns
        -------
        int
            The continous mode state.
        """        
        cmd = lambda c: b":INIT:GATE" + str(c).encode() 
        buffer = 32
        set_str = state.encode()

        if chn == 0:
            self._sendSCPI(cmd("") + b" " + set_str +  b"\n")
        elif chn == 3:
            self._sendSCPI(cmd(1) + b" " + set_str +  b"\n")
            self._sendSCPI(cmd(2) + b" " + set_str +  b"\n")
            
        elif chn <3:
            self._sendSCPI(cmd(chn) + str(chn).encode + b" " + set_str+ "\n")

        return self.getContMode(chn)

    def getTrigLevel(self) -> float:
        """getTrigLevel Get the trigger level in Volts

        Returns
        -------
        float
            The trigger level in volts
        """        
        cmd = lambda c: b":ARM:TRIG:LEV" + str(c).encode()
        buffer = 32
        response = self._sendrecvSCPI(cmd("") + b"?" +  b"\n",buffer)
        
        return float(response.decode())

    def setTrigLevel(self,value:float or int or str):
        """setTrigLevel Sets the Trigger Level in Volts

        Parameters
        ----------
        value : float or int or str
            The value of the trigger level. Either in Volts or "min"/"max"

        Returns
        -------
        float
            The new trigger level in volts

        """        
        if type(value) == str:
            if value.lower() == "min":
                value_str = b"MIN"
            elif value.lower() == "max":
                value_str = b"MAX"
            else:
                raise ValueError("If value is a string it must be either 'MIN' or 'MAX'")
        else:
            value_str = str(value).encode()

        cmd = lambda c: b":ARM:TRIG:LEV" + str(c).encode()
        self._sendSCPI(cmd("") + b" " + value_str +  b"\n")

        return self.getTrigLevel()

    def getTrigImp(self) -> str:
        """getTrigImp Get the trigger impedance.

        Returns
        -------
        str
            The Trigger Impedence. Either: High or Low.
        """        
        cmd = lambda c: b":ARM:TRIG:IMP" + str(c).encode()
        buffer = 32
        response = self._sendrecvSCPI(cmd("") + b"?" +  b"\n",buffer)
        
        return str(response.decode())

    def setTrigImp(self,level:str) -> str:
        """setTrigImp Set the trigger impedence

        Parameters
        ----------
        level : str
            The Trigger Impedence. Either: High or Low.

        Returns
        -------
        str
            The Trigger Impedence. Either: High or Low.
        """
        
        if level.lower() == "low":
            value_str = b"LOW"
        elif level.lower() == "high":
            value_str = b"HIGH"
        
        cmd = lambda c: b":ARM:TRIG:IMP" + str(c).encode()
        self._sendSCPI(cmd("") + b" " + value_str +  b"\n")

        return self.getTrigImp()