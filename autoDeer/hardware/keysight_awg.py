"""
This is selection of functions and comands for the control of the keysight M8120 AWG through a TCP Socket
"""
import os.path as path
import re
import socket
from typing import Optional,Union
import numpy as np
from io import StringIO
import logging
import matplotlib.pyplot as plt

hw_log = logging.getLogger('hardware.AWG')

def dacDecode(data):
    data_s = data.decode()
    raw_values = np.genfromtxt(StringIO(data_s), delimiter = ",")
    num_values = len(raw_values)
    markers = np.zeros(2,num_values)
    decoded_data = np.zeros(num_values)
    for i,val in enumerate(raw_values):
        markers[0,i] = val & 1
        markers[1,i] = (val & 2) >> 1
        decoded_data[i] = val // 16
    return decoded_data


def dacEncode(data,markers):
    num_values = len(data)
    num_markers = np.shape(markers)[0]
    new_data = np.zeros(num_values,dtype=np.int16)
    for i in range(0,num_values):
        tmp_data = data[i] * 16
        tmp_data = tmp_data + int(markers[0,i])
        tmp_data = tmp_data + int(markers[1,i] * 2)
        new_data[i] = tmp_data
    return new_data
    
class sequence:

    def __init__(self,awg) -> None:
        self.awg = awg
        pass

class waveform:

    def __init__(self,awg,shape,markers) -> None:
        self.awg = awg
        self.shape = shape
        self.mk1 = markers[0,:]
        self.mk2 = markers[1,:]

        pass

    def _awg_check(self):
        sampling_freq = 12e9
        grad:int = 64
        min_length:int = 5*grad
        max = 2**11

        # Check_amplitude
        if np.max(np.abs(self.shape)) >= max:
            raise WaveformError(f"Waveform exceeds maximum possible value. Limited to +- {max}")

        # Check Gradularity
        if np.shape(self.shape)[0]%grad != 0:
            raise WaveformError(f"Waveform gradularity must equal {grad}")

        # Check Minimum Length
        if np.shape(self.shape)[0] < min_length:
            raise WaveformError(f"Waveform must be longer than {min_length}")

        return True

    def set_amp(self,val:Union[float,str]):
        max = 2**11 - 1

        if type(val) == str:
            if val.lower() == 'max':
                val = 1
            elif val.lower == 'min':
                val = 0
        
        cur_max = np.max(abs(self.shape))
        self.shape = self.shape / cur_max
        self.shape = self.shape * max
        self.shape = self.shape.astype(np.int16)
 
    def DACencode(self):
        encodedData = dacEncode(self.shape,np.vstack([self.mk1,self.mk2]))
        return encodedData

    def plot(self):
        fig, ax = plt.subplots()
        ax.plot(self.shape)
        ax.set_xlabel('Sample')
        ax.set_ylabel('DAC')
        return fig


class interface:

    def __init__(self) -> None:
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

            # self.instr_socket.sendall(b"*RST;*OPC?\n")
            # result = self.instr_socket.recv(8)

            self.instr_socket.sendall(b"*IDN?\n")
            idn_result = self.instr_socket.recv(255)
            print(f"Identification String: {idn_result}")

        except socket.error:
            print("Instrument initalisation failed")
            self.instr_socket.close()
        pass

    def _sendSCPI(self,cmd,cmd_check:bool=True,error_check:bool = True) -> None:

        if type(cmd) is str:
            cmd_b = cmd.encode()
            cmd_s = cmd
        elif type(cmd) is bytes:
            cmd_b = cmd
            cmd_s = cmd.decode()

        if cmd_check == True:
            self._check_cmd(cmd_s)
        
        self.instr_socket.sendall(cmd_b)
        hw_log.info(f'SCPI_Cmd:{cmd_s}')

        if error_check:
            self._raise_errors()
        
        pass

    def _recvSCPI(self,buffer,msg_check:bool = True,decode:bool = True):
        data_b = self.instr_socket.recv(buffer)

        if decode:
            data = data_b.decode()
            hw_log.debug(f'SCPI_recv:{data}')

            if msg_check:
                self._check_msg(data)
            
            return  data
        else:
            return data_b

    def _sendrecvSCPI(self,cmd:str,buffer:int,error_check:bool=True,cmd_check:bool = True):

        self._sendSCPI(cmd,cmd_check=cmd_check,error_check=False)
        data = self._recvSCPI(buffer)

        if error_check:
            self._raise_errors()
        
        return  data

    def _get_errors(self):

        errors = []
        check = True
        while check:
            error = self._sendrecvSCPI(b":SYST:ERR? \n",512,error_check=False)
            if error[0] == '0':
                check = False
            else:
                errors.append(error)
        return errors

    def _check_errors(self):

        errors = self._get_errors()
        if errors == []:
            return False
        else:
            return True

    def _raise_errors(self):
        errors = self._get_errors()

        if errors == []:
            pass
        else:
            msg = ""
            for e in errors:
                msg += e
            raise RuntimeError(msg)

    def _check_cmd(self,cmd:str or bytes):

        # Check the line termintation is included

        if type(cmd) == bytes:
            cmd_s = cmd.decode()
        else:
            cmd_s = cmd

        check_string = "\n\Z"
        match = re.search(check_string,cmd_s)
        
        if match == None:
            raise ValueError("The given SCPI AWG cmd is not terminated (\\n) \n cmd:"+cmd_s)

    def _check_msg(self,msg:str or bytes):

        # Check the line termintation is included

        if type(msg) == bytes:
            msg_s = msg.decode()
        else:
            msg_s = msg

        check_string = "\n\Z"
        match = re.search(check_string,msg_s)
        
        if match == None:
            raise ValueError("The recieved TCP msg is not terminated: Please increase the buffer size")

    def getSystemSettings(self,savefile:Optional[str]=None):
        """SystemSettings Gets the current AWG settings. This is either returns as a string or saves to file.

        Parameters
        ----------
        savefile : str or None, optional
            The file path of where to save the system settings. If None is give returns a string, by default None

        Returns
        -------
        Optional[str]
            Returns the system settings as one long string.
        """        
        data = self._sendrecvSCPI(b':SYST:SET? \n',20000,error_check=True)

        # This removes any #xxxx that occurs at the begining. This can't be loaded by the spectrometer.
        re_mask = "^#[0-9]*"
        data = re.sub(re_mask,"",data)
        
        if savefile == None:
            return data
        else:
            with open(savefile, 'w') as f:
                f.write(data)
            pass

    def setSystemSettings(self,settings:Optional[Union[bytes,str]] = None,save_path:Optional[str] = None)->None:
        """SystemSettings Sends and sets the AWG settings. This can either take a string of settings or a save file

        Parameters
        ----------
        settings : Optional[Union[bytes,str]]
            This can be either a bytes string or a normal string of the settings.
        save_path : Optional[str]
            The file path to a save file containg the settings as a single line string.
        """

        if (not settings) and (not save_path):
            raise ValueError('Only one of settings or save_path can be given.')
        elif settings != None:
            if type(settings) == str:
                data = settings.encode()
            else:
                data = settings
        elif save_path != None:
            if path.exists(save_path):
                with open(save_path,'r') as f:
                    data = f.read()
                data = data.encode()
        else:
            raise ValueError('One of settings or save_path can be given.')

        msg = b':SYST:SET'
        msg += data
        self._sendSCPI(msg,cmd_check=True)
        self._raise_errors()


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


class WaveformError(Exception):
    pass

class AWGError(Exception):
    pass

class SequenceError(Exception):
    pass