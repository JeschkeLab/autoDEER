"""
This is selection of functions and comands for the control of the keysight M8120 AWG through a TCP Socket
"""
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


    def Abort(self,chn:int):
        """Abort Stop signal generation on either one or more channels. If channels are coupled bith channels are stopped.

        Parameters
        ----------
        chn : int
            The channel number to be stopped. 1 or 2 stops the respective channel and 3 stops both
        """        
        if chn == 3:
            self._sendSCPI(b":ABOR1\n")
            self._sendSCPI(b":ABOR2\n")
        elif chn <3:
            self.__sendSCPI(b":ABOR"+ str(chn).encode + "\n")
    
