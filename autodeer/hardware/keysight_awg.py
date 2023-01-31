"""
This is selection of functions and comands for the control of the keysight
M8120 AWG through a TCP Socket
"""
import os.path as path
import re
import socket
from typing import Optional, Union
import numpy as np
from io import StringIO
import logging
import matplotlib.pyplot as plt
hw_log = logging.getLogger('hardware.AWG')

# ============================================================================


def dacDecode(data):
    data_s = data.decode()
    raw_values = np.genfromtxt(StringIO(data_s), delimiter=",")
    num_values = len(raw_values)
    markers = np.zeros(2, num_values)
    decoded_data = np.zeros(num_values)
    for i, val in enumerate(raw_values):
        markers[0, i] = val & 1
        markers[1, i] = (val & 2) >> 1
        decoded_data[i] = val // 16
    return decoded_data

# ============================================================================


def dacEncode(data, markers):
    num_values = len(data)
    # num_markers = np.shape(markers)[0]
    new_data = np.zeros(num_values, dtype=np.int16)
    for i in range(0, num_values):
        tmp_data = data[i] * 16
        tmp_data = tmp_data + int(markers[0, i])
        tmp_data = tmp_data + int(markers[1, i] * 2)
        new_data[i] = tmp_data
    return new_data

# ============================================================================


class SequenceTable:

    def __init__(self, awg) -> None:
        self.awg = awg
        pass
    
    def read_from_str(self, str):
        re_mask = "[0-9]+,[0-9]+,[0-9]+,[0-9]+,[0-9]+,[0-9]+"
        entries = re.findall(re_mask, str)
        self.table = []
        for entry in entries:
            new_seq = Sequence(self.awg)
            new_seq.read_from_str(entry)
            self.table.append(new_seq)
        return self.table

    def write_to_str(self):
        string = ""
        for seq in self.table:
            string += seq.write_to_str()
            string += ","
        string = string[:-1]

        return string
    

# ============================================================================


class Sequence:

    def __init__(self, awg) -> None:
        self.awg = awg
        pass

    def read_from_str(self, str):
        re_mask = "[0-9]+"
        parameters = re.findall(re_mask, str)
        self.params = []
        for param in parameters:
            self.params.append(int(param))

    def write_to_str(self):
        string = ""
        for param in self.params:
            string += (str(param) + ",")
        string = string[:-1]
        return string
        
    def decode(self, mode):
        self.entry_dict = {}
        for i, param in enumerate(self.params):
            param = int(param)

            if i == 0:
                data_cmd = (param & 2**31) // 2**31
                self.entry_dict['data_cmd'] = data_cmd
                self.entry_dict['end_marker_seq'] = (param & 2**30)//2**30
                self.entry_dict['end_marker_sce'] = (param & 2**29)//2**29
                self.entry_dict['init_marker_seq'] = (param & 2**28)//2**28
                self.entry_dict['marker_enable'] = (param & 2**24)//2**24
                # 0b111100000000000000000000
                self.entry_dict['seq_adv_mode'] = (param & 15728640) // 2**20
                # 0b11110000000000000000
                self.entry_dict['seg_adv_mode'] = (param & 983040)//2**16
                self.entry_dict['amp_tab_init'] = (param & 2**15)//2**15
                self.entry_dict['amp_tab_incr'] = (param & 2**14)//2**14
                self.entry_dict['freq_tab_init'] = (param & 2**13)//2**13
                self.entry_dict['freq_tab_incr'] = (param & 2**12)//2**12

                if data_cmd == 0:
                    # We can identify that this is a Data entry
                    self.entry_type = 'data'
                    self.entry_dict['entry_type'] = self.entry_type
                else:
                    # We still need the cmd code in word 3 to identify if it is 
                    # idle delay or config
                    self.entry_type = None
                
            elif i == 1:
                # The sequence loop count. This should be zero in all situation
                #  except the first entry of a sequence
                self.entry_dict['seq_loop_cnt'] = param

            elif i == 2:

                if self.entry_type == 'data':
                    self.entry_dict['seg_loop_cnt'] = param
        
                else:
                    if (param & 65535) == 0:  # 0b1111111111111111
                        self.entry_type = 'idle delay'
                        self.entry_dict['entry_type']
                        self.entry_dict['cmd_code'] = 0

                    elif (param & 65535) == 1:  # 0b1111111111111111
                        self.entry_type = 'config'
                        self.entry_dict['entry_type'] = self.entry_type
                        self.entry_dict['cmd_code'] = 1
                        self.entry_dict['action_table_id'] = \
                            (param & 4294901760) // 2**16

            elif i == 3:

                if self.entry_type == 'idle delay':
                    if mode == 'speed':
                        self.entry_dict['idle_sample'] = (param & 4095)
                    else:
                        raise RuntimeError(
                            "Only Speed mode has been currently implemented "
                            "for the idle delay cmd.")

                elif (self.entry_type == 'data') or \
                     (self.entry_type == 'config'):
                    self.entry_dict['seg_id'] = (param & 0b1111111111111111111)

                else:
                    raise RuntimeError(
                        f"entry_type must be one of"
                        f"['data','idle delay','config']. Currently it is"
                        f" {self.entry_type}")     
                
            elif i == 4:

                if self.entry_type == 'idle delay':
                    
                    self.entry_dict['idle_delay'] = param 

                elif (self.entry_type == 'data') or \
                     (self.entry_type == 'config'):
                    self.entry_dict['seg_start_offset'] = param

                else:
                    raise RuntimeError(
                        f"entry_type must be one of ['data','idle delay',"
                        f"'config']. Currently it is {self.entry_type}")   
                    
            elif i == 5:

                if self.entry_type == 'idle delay':
                    param = 0 
                elif (self.entry_type == 'data') or \
                     (self.entry_type == 'config'):
                    self.entry_dict['seg_end_offset'] = param

                else:
                    raise RuntimeError(
                        f"entry_type must be one of ['data','idle delay',"
                        f"'config']. Currently it is {self.entry_type}")    

    def _binary_builder(self, dict, element, bit_pos):
        if element in dict:
            return dict[element] * 2**bit_pos
        else:
            return 0

    def encode(self, dict):
        # mode = 'speed'
        # Building the control word
        control_word = 0
        control_word += self._binary_builder(dict, 'data_cmd', 31)
        control_word += self._binary_builder(dict, 'end_marker_seq', 30)
        control_word += self._binary_builder(dict, 'end_marker_sce', 29)
        control_word += self._binary_builder(dict, 'init_marker_seq', 28)
        control_word += self._binary_builder(dict, 'marker_enable', 24)
        control_word += self._binary_builder(dict, 'seq_adv_mode', 20)
        control_word += self._binary_builder(dict, 'seg_adv_mode', 16)
        control_word += self._binary_builder(dict, 'amp_tab_init', 15)
        control_word += self._binary_builder(dict, 'amp_tab_incr', 14)
        control_word += self._binary_builder(dict, 'freq_tab_init', 13)
        control_word += self._binary_builder(dict, 'freq_tab_incr', 12)
        # print(control_word)

        # Building the Segment loop word
        seq_loop_word = dict['seq_loop_cnt']

        if dict['data_cmd'] == 0:
            word2 = dict['seg_loop_cnt']
            word3 = dict['seg_id']
            word4 = dict['seg_start_offset']
            word5 = dict['seg_end_offset']
        elif dict['data_cmd'] == 1:
            word2 = dict['cmd_code']
            if dict['cmd_code'] == 0:
                word3 = dict['idle_sample']
                word4 = dict['idle_delay']
                word5 = 0

            elif dict['cmd_code'] == 1:
                word2 += self._binary_builder(dict, 'action_table_id', 16)
                word3 = dict['seg_id']
                word4 = dict['seg_start_offset']
                word5 = dict['seg_end_offset']
        entry = str(control_word) + ',' + str(seq_loop_word) + ',' + \
            str(word2) + ',' + str(word3) + ',' + str(word4) + ',' + str(word5)
        self.params = [control_word, seq_loop_word, word2, word3, word4, word5]
        return entry

# =============================================================================


class Waveform:

    def __init__(self, awg, shape, markers=None) -> None:
        self.awg = awg

        if np.ndim(shape) == 2:
            self.complex = True
        
        self.shape = shape
        self.num_points = np.shape(self.shape)
        if type(markers) is None:
            self.mk1 = np.ones(self.num_points, dtype=np.int16)
            self.mk2 = np.ones(self.num_points, dtype=np.int16)
        else:
            self.mk1 = markers[0, :]
            self.mk2 = markers[1, :]

        pass

    def _awg_check(self):
        # sampling_freq = 12e9
        grad: int = 64
        min_length: int = 5*grad
        shape_max = 2**11
        encode_max = 2**15

        # Check amplitude of waveform
        if np.max(np.abs(self.shape)) > shape_max:
            raise WaveformError(
                f"Waveform exceeds maximum possible value. Limited to"
                f" +- {shape_max}")

        # Check Gradularity
        if self.complex is True:
            n = np.shape(self.shape)[1]
        else:
            n = np.shape(self.shape)[0]

        if n % grad != 0:
            raise WaveformError(f"Waveform gradularity must equal {grad}")

        # Check Minimum Length
        if n < min_length:
            raise WaveformError(f"Waveform must be longer than {min_length}")

        # Check encoded data
        if hasattr(self, 'encodedData'):
            if self.complex is True:
                en = np.shape(self.encodedData)[1]
            else:
                en = np.shape(self.encodedData)[0]
                
            if np.max(np.abs(self.encodedData)) > encode_max:
                raise WaveformError(
                    f"Waveform exceeds maximum possible value."
                    f"Limited to +- {encode_max}")
            # Check Gradularity
            if en % grad != 0:
                raise WaveformError(f"Waveform gradularity must equal {grad}")

            # Check Minimum Length
            if en < min_length:
                raise WaveformError(
                    f"Waveform must be longer than {min_length}")

        return True

    def set_amp(self, val: Union[float, str]):
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
        encodedData_real = dacEncode(
            self.shape[0, :], np.vstack([self.mk1, self.mk2]))
        encodedData_imag = dacEncode(
            self.shape[0, :], np.vstack([self.mk1, self.mk2]))
        self.encodedData = np.vstack([encodedData_real, encodedData_imag])
        return self.encodedData

    def write_to_string(self):
        self.DACencode()

        if self.complex is True:
            string1 = ""
            for num in self.encodedData[0, :]:
                string1 += str(int(num))
                string1 += ","
            string1 = string1[:-1]    

            string2 = ""
            for num in self.encodedData[1, :]:
                string2 += str(int(num))
                string2 += ","
            string2 = string2[:-1] 

            return [string1, string2]
        else:
            string = ""
            for num in self.encodedData:
                string += str(int(num))
                string += ","
            string = string[:-1]
            return string

    def plot(self):
        fig, ax = plt.subplots()
        if self.complex is True:
            ax.plot(self.shape[0, :])
            ax.set_title('Channel 1')
        else:
            ax.plot(self.shape)
        ax.set_xlabel('Sample')
        ax.set_ylabel('DAC')
        return fig

    def pulse_blanking(self, preamb=1800, postamb=0):
        # Blank Channels
        if self.complex is True:
            prestack = np.zeros((2, preamb), dtype=np.int16)
            poststack = np.zeros((2, postamb), dtype=np.int16)
        else:
            prestack = np.zeros(preamb, dtype=np.int16)
            poststack = np.zeros(postamb, dtype=np.int16)

        self.shape = np.hstack([prestack, self.shape, poststack])

        # Blank Markers
        if self.complex is True:
            prestack = np.ones((2, preamb), dtype=np.int16)
            poststack = np.ones((2, postamb), dtype=np.int16)
        else:
            prestack = np.ones(preamb, dtype=np.int16)
            poststack = np.ones(postamb, dtype=np.int16)
        self.mk1 = np.hstack([prestack, self.mk1, poststack])
        self.mk2 = np.hstack([prestack, self.mk2, poststack])
        
        self.num_points = self.shape.shape[-1]

    def enforce_gradualrity(self):
        if self.complex is True:
            wave_length = np.shape(self.shape)[1]
            if (wave_length % 64) != 0:
                extra_zeros = np.zeros((2, 64-wave_length % 64))
                self.shape = np.hstack([self.shape, extra_zeros])

                extra_ones = np.ones((2, 64-wave_length % 64))
                self.mk1 = np.hstack([self.mk1, extra_ones])
                self.mk2 = np.hstack([self.mk2, extra_ones])
        else:
            wave_length = np.shape(self.shape)[0]
            if (wave_length % 64) != 0:
                extra_zeros = np.zeros(64 - wave_length % 64)
                self.shape = np.hstack([self.shape, extra_zeros])

                extra_ones = np.ones(64-wave_length % 64)
                self.mk1 = np.hstack([self.mk1, extra_ones])
                self.mk2 = np.hstack([self.mk2, extra_ones])

        self.num_points = self.shape.shape[-1]

    def define_new_waveform(self, ID: int, ch: int = 3):

        if self.complex is True:
            # print(f"Complex Waveform so ignoring channel. Real-> Ch1, 
            # Imag -> Ch2")
            self.awg.defineWaveform(ID, self.num_points, 0, ch=3)
        else:
            self.awg.defineWaveform(ID, self.num_points, 0, ch=ch)

    def send_waveform(self, ID: int, ch: int = 3):
        
        string = self.write_to_string()
        if self.complex is True:
            # print(f"Complex Waveform so ignoring channel. Real-> Ch1, 
            # Imag -> Ch2")
            
            self.awg.setWaveform(string[0], ID, 1, 0)
            self.awg.setWaveform(string[1], ID, 2, 0)

        else:
            self.awg.setWaveform(string, ID, ch, 0)


class Interface:

    def __init__(self) -> None:
        pass

    def open_con(self, host, port=5025) -> None:
        try:
            self.instr_socket = socket.socket(
                socket.AF_INET, socket.SOCK_STREAM)
        except socket.error as msg:
            print('Failed to create socket. Error code: ' + str(msg[0]) +
                  ' , Error message : ' + msg[1])
            self.instr_socket.close()
        
        print("Socket Created")

        self.wait_time = 0.2  # A delay time variable in seconds
        
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
        self.instr_socket.connect((self.ip, self.port))

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

    def _sendSCPI(
            self, cmd, cmd_check: bool = True,
            error_check: bool = True) -> None:

        if type(cmd) is str:
            cmd_b = cmd.encode()
            cmd_s = cmd
        elif type(cmd) is bytes:
            cmd_b = cmd
            cmd_s = cmd.decode()
        else:
            raise ValueError(
                f"Cmd must be either of type str or bytes. The type "
                f"is:{type(cmd)}")

        if cmd_check is True:
            self._check_cmd(cmd_s)
        
        self.instr_socket.sendall(cmd_b)
        hw_log.info(f'SCPI_Cmd:{cmd_s}')

        if error_check:
            self._raise_errors()
        
        pass

    def _recvSCPI(self, buffer, msg_check: bool = True, decode: bool = True):
        data_b = self.instr_socket.recv(buffer)

        if decode:
            data = data_b.decode()
            hw_log.debug(f'SCPI_recv:{data}')

            if msg_check:
                self._check_msg(data)
            return data

        else:
            return data_b

    def _sendrecvSCPI(
            self, cmd: str, buffer: int, decode: int = True,
            error_check: bool = True, cmd_check: bool = True):

        self._sendSCPI(cmd, cmd_check=cmd_check, error_check=False)
        data = self._recvSCPI(buffer, decode)

        if error_check:
            self._raise_errors()
        
        return data

    def _get_errors(self):

        errors = []
        check = True
        while check:
            error = self._sendrecvSCPI(
                b":SYST:ERR? \n", 512, error_check=False)
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

    def _check_cmd(self, cmd: str or bytes):

        # Check the line termintation is included

        if type(cmd) == bytes:
            cmd_s = cmd.decode()
        else:
            cmd_s = cmd

        check_string = r"\n\Z"
        match = re.search(check_string, cmd_s)
        
        if match is None:
            raise ValueError(
                "The given SCPI AWG cmd is not terminated (\\n) \n cmd:"+cmd_s)

    def _check_msg(self, msg: str or bytes):

        # Check the line termintation is included

        if type(msg) == bytes:
            msg_s = msg.decode()
        else:
            msg_s = msg

        check_string = r"\n\Z"
        match = re.search(check_string, msg_s)
        
        if match is None:
            raise ValueError(
                "The recieved TCP msg is not terminated: Please increase the"
                " buffer size")

    def getSystemSettings(self, savefile: Optional[str] = None):
        """SystemSettings Gets the current AWG settings. This is either 
        returns as a string or saves to file.

        Parameters
        ----------
        savefile : str or None, optional
            The file path of where to save the system settings. If None is 
            give returns a string, by default None

        Returns
        -------
        Optional[str]
            Returns the system settings as one long string.
        """        
        data = self._sendrecvSCPI(b':SYST:SET? \n', 20000, error_check=True)

        # This removes any #xxxx that occurs at the begining. 
        # This can't be loaded by the spectrometer.
        
        re_mask = "^#[0-9]*"
        data = re.sub(re_mask, "", data)
        
        if savefile is None:
            return data
        else:
            with open(savefile, 'w') as f:
                f.write(data)
            pass

    def setSystemSettings(
            self, settings: Optional[Union[bytes, str]] = None,
            save_path: Optional[str] = None) -> None:
        """SystemSettings Sends and sets the AWG settings. This can either 
        take a string of settings or a save file

        Parameters
        ----------
        settings : Optional[Union[bytes,str]]
            This can be either a bytes string or a normal string of the 
            settings.
        save_path : Optional[str]
            The file path to a save file containg the settings as a single 
            line string.
        """

        if (not settings) and (not save_path):
            raise ValueError('Only one of settings or save_path can be given.')
        elif settings is not None:
            if type(settings) == str:
                data = settings.encode()
            else:
                data = settings
        elif save_path is not None:
            if path.exists(save_path):
                with open(save_path, 'r') as f:
                    data = f.read()
                data = data.encode()
        else:
            raise ValueError('One of settings or save_path can be given.')

        msg = b':SYST:SET'
        msg += data
        self._sendSCPI(msg, cmd_check=True)
        self._raise_errors()

    def Abort(self, chn: int = 0) -> None:
        """Abort Stop signal generation on either one or more channels. If 
        channels are coupled bith channels are stopped.

        Parameters
        ----------
        chn : int
            The channel number. 1 or 2 for the respective channel, 3 for both 
            and 0 for default/coupled, by default =0
        """        
        if chn == 3:
            self._sendSCPI(b":ABOR1\n")
            self._sendSCPI(b":ABOR2\n")
        elif chn == 0:
            self._sendSCPI(b":ABOR\n")
        elif chn < 3:
            self.__sendSCPI(b":ABOR" + str(chn).encode + "\n")

    def StartSigGen(self, chn: int = 0) -> None:
        """StartSigGen Start signal generation on a channel. Of channels are 
        coupled, both channels are started

        Parameters
        ----------
        chn : int
            The channel number. 1 or 2 for the respective channel, 3 for both 
            and 0 for default/coupled, by default =0
        """
        if chn == 3:
            self._sendSCPI(b":INIT:IMM1\n")
            self._sendSCPI(b":INIT:IMM2\n")
        elif chn == 0:
            self._sendSCPI(b":INIT:IMM\n")
        elif chn < 3:
            self.__sendSCPI(b":INIT:IMM1" + str(chn).encode + "\n")
        
        pass
    
    def getTrigDelay(self, chn: int = 0) -> float:
        """getTrigDelay Get the trigger to data out delay.

        Parameters
        ----------
        chn : int, optional
            The channel number. 1 or 2 for the respective channel and 0 for 
            default/coupled, by default =0

        Returns
        -------
        float
            Read the trigger delay in seconds
        """        
        cmd = b":ARM:DEL"
        buffer = 32
        if chn == 0:
            response = self._sendrecvSCPI(cmd + b"?" + b"\n", buffer)
        elif chn < 3:
            response = self._sendrecvSCPI(
                cmd + str(chn).encode + b"?" + "\n", buffer)

        return float(response.decode())

    def setTrigDelay(self, delay: float, chn: int = 0) -> float:
        """setTrigDelay Set the trigger to data out delay.

        Parameters
        ----------
        delay : float
            The trigger delay in seconds
        chn : int, optional
            The channel number. 1 or 2 for the respective channel, 3 for both 
            and 0 for default/coupled, by default =0

        Returns
        -------
        float
            Read the trigger delay in seconds
        """

        cmd = b":ARM:DEL"
        set_str = str(delay).encode()

        if chn == 3:
            self._sendSCPI(cmd + b"1 " + set_str + b"\n")
            self._sendSCPI(cmd + b"2 " + set_str + b"\n")

        if chn == 0:
            self._sendSCPI(cmd + b" " + set_str + b"\n")
        elif chn < 3:
            self._sendSCPI(
                cmd + str(chn).encode + b" " + set_str + "\n")
        
        return self.getTrigDelay(chn)

    def getArmMode(self, chn: int = 0) -> str:
        """getArmMode Get the arming mode

        Parameters
        ----------
        chn : int, optional
            The channel number. 1 or 2 for the respective channel and 0 for 
            default/coupled, by default =0

        Returns
        -------
        str
            Read the arming mode
        """        
        cmd = lambda c: b":INIT:CONT" + str(c).encode() + b":ENABI"
        buffer = 32
        if chn == 0:
            response = self._sendrecvSCPI(cmd("") + b"?" + b"\n", buffer)
        elif chn < 3:
            response = self._sendrecvSCPI(
                cmd(chn) + str(chn).encode + b"?" + "\n", buffer)

        return response.decode()

    def setArmMode(self, mode: str, chn: int = 0) -> str:
        """setArmMode Set the arming mode

        Parameters
        ----------
        mode : str
            The arming mode. "SELF" or "ARM"ed
        chn : int, optional
            The channel number. 1 or 2 for the respective channel, 3 for both 
            and 0 for default/coupled, by default =0

        Returns
        -------
        str
            Returns the state of the arming mode
        """        
        
        cmd = lambda c: b":INIT:CONT" + str(c).encode() + b":ENABI"
        set_str = mode.encode()

        if chn == 0:
            self._sendSCPI(cmd("") + b" " + set_str + b"\n")
        elif chn == 3:
            self._sendSCPI(cmd(1) + b" " + set_str + b"\n")
            self._sendSCPI(cmd(2) + b" " + set_str + b"\n")
            
        elif chn < 3:
            self._sendSCPI(cmd(chn) + str(chn).encode + b" " + set_str + "\n")

        return self.getArmMode(chn)

    def getContMode(self, chn: int = 0) -> int:
        """getContMode Get the continous mode.

        Parameters
        ----------
        chn : int, optional
            The channel number. 1 or 2 for the respective channel and 0 for 
            default/coupled, by default =0

        Returns
        -------
        int
            The continous mode state.
        """        

        cmd = lambda c: b":INIT:CONT" + str(c).encode()
        buffer = 32
        if chn == 0:
            response = self._sendrecvSCPI(cmd("") + b"?" + b"\n", buffer)
        elif chn < 3:
            response = self._sendrecvSCPI(
                cmd(chn) + str(chn).encode + b"?" + "\n", buffer)

        return int(response.decode())

    def setContMode(self, state: int, chn: int = 0) -> int:
        """getContMode Set the continous mode. This must be used in conjuction
         with "setGateMode".

        Parameters
        ----------
        state : int
            The state of the continous mode. 0 for off, 1 for on.
        chn : int, optional
            The channel number. 1 or 2 for the respective channel, 3 for both 
            and 0 for default/coupled, by default =0

        Returns
        -------
        int
            The continous mode state.
        """        
        cmd = lambda c: b":INIT:CONT" + str(c).encode() 
        set_str = state.encode()

        if chn == 0:
            self._sendSCPI(cmd("") + b" " + set_str + b"\n")
        elif chn == 3:
            self._sendSCPI(cmd(1) + b" " + set_str + b"\n")
            self._sendSCPI(cmd(2) + b" " + set_str + b"\n")

        elif chn < 3:
            self._sendSCPI(cmd(chn) + str(chn).encode + b" " + set_str + "\n")

        return self.getContMode(chn)

    def getGateMode(self, chn: int = 0) -> int:
        """getGateMode Get the continous mode.

        Parameters
        ----------
        chn : int, optional
            The channel number. 1 or 2 for the respective channel and 0 for 
            default/coupled, by default =0

        Returns
        -------
        int
            The Gate mode.
        """        

        cmd = lambda c: b":INIT:GATE" + str(c).encode()
        buffer = 32
        if chn == 0:
            response = self._sendrecvSCPI(cmd("") + b"?" + b"\n", buffer)
        elif chn < 3:
            response = self._sendrecvSCPI(
                cmd(chn) + str(chn).encode + b"?" + "\n", buffer)

        return int(response.decode())

    def setGateMode(self, state: int, chn: int = 0) -> int:
        """getContMode Set the continous mode. This must be used in conjuction 
        with "setContMode".

        Parameters
        ----------
        state : int
            The state of the gate mode. 0 for off, 1 for on.
        chn : int, optional
            The channel number. 1 or 2 for the respective channel, 3 for both 
            and 0 for default/coupled, by default =0

        Returns
        -------
        int
            The continous mode state.
        """        
        cmd = lambda c: b":INIT:GATE" + str(c).encode() 
        set_str = state.encode()

        if chn == 0:
            self._sendSCPI(cmd("") + b" " + set_str + b"\n")
        elif chn == 3:
            self._sendSCPI(cmd(1) + b" " + set_str + b"\n")
            self._sendSCPI(cmd(2) + b" " + set_str + b"\n")
            
        elif chn < 3:
            self._sendSCPI(cmd(chn) + str(chn).encode + b" " + set_str + "\n")

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
        response = self._sendrecvSCPI(cmd("") + b"?" + b"\n", buffer)
        
        return float(response.decode())

    def setTrigLevel(self, value: float or int or str):
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
                raise ValueError(
                    "If value is a string it must be either 'MIN' or 'MAX'")
        else:
            value_str = str(value).encode()

        cmd = lambda c: b":ARM:TRIG:LEV" + str(c).encode()
        self._sendSCPI(cmd("") + b" " + value_str + b"\n")

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
        response = self._sendrecvSCPI(cmd("") + b"?" + b"\n", buffer)
        
        return str(response.decode())

    def setTrigImp(self, level: str) -> str:
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
        self._sendSCPI(cmd("") + b" " + value_str + b"\n")

        return self.getTrigImp()

    def getMemoryStatus(self, ch: int = 3) -> np.ndarray:

        def print_mem(mem, ch):
            re_mask = "[0-9]+"
            match = re.findall(re_mask, mem)
            print(
                f"Channel {ch}: Bytes Avaliable = {int(match[0])},Bytes in "
                f"use = {int(match[1])},Continguous Bytes Avaliable = "
                f"{int(match[3])}")

            return [int(match[0]), int(match[1]), int(match[2])]

        if ch == 3:
            mem1 = self._sendrecvSCPI(b"TRAC1:FREE? \n", 128, decode=True)
            mem2 = self._sendrecvSCPI(b"TRAC2:FREE? \n", 128, decode=True)
            mem1 = print_mem(mem1, 1)
            mem2 = print_mem(mem2, 2)
            mem = np.vstack(mem1, mem2)

        elif (ch == 1) or (ch == 3):
            cmd = "TRAC{}:FREE? \n".format(int(ch))
            mem1 = self._sendrecvSCPI(cmd.encode(), 128, decode=True)
            mem = print_mem(mem1, ch)

        else:
            raise ValueError("Channel must be one of 1,2,3(both).")

        return mem

    def getWaveformCatalog(self, ch: int = 3) -> np.ndarray:

        def extractCat(cat_string):
            cat_mask = "[0-9]+,[0-9]+"
            waveforms = re.findall(cat_mask, cat_string)
            num_waveforms = len(waveforms)
            catalog = np.zeros((num_waveforms, 2), dtype=np.int16)
            for i, wave in enumerate(waveforms):
                iD, num_points = re.findall("[0-9]+", wave)
                catalog[i, :] = [int(iD), int(num_points)]
            return catalog

        if ch == 3:
            
            cat1 = self._sendrecvSCPI(b"TRAC1:CAT? \n", 128)
            cat2 = self._sendrecvSCPI(b"TRAC2:CAT? \n", 128)
            return [extractCat(cat1), extractCat(cat2)]

        elif (ch == 1) or (ch == 2):
            
            cmd = "TRAC{}:CAT? \n".format(int(ch))
            cat1 = self._sendrecvSCPI(cmd.encode(), 128)
            return extractCat(cat1)

        else:
            raise ValueError("Channel must be one of 1,2,3(both).")

    def deleteWaveform(self, ch: int, ID: int) -> None:

        if ch == 3:
            cmd = "TRAC1:DEL {} \n".format(ID)
            self._sendSCPI(cmd.encode)
            cmd = "TRAC2:DEL {} \n".format(ID)
            self._sendSCPI(cmd.encode)
        
        elif (ch == 1) or (ch == 2):
            cmd = "TRAC{}:DEL {} \n".format(ch, ID)
            self._sendSCPI(cmd.encode())
        
        else:
            raise ValueError("Channel must be one of 1,2,3(both).")

        pass

    def deleteAllWaveforms(self, ch: int) -> None:

        if ch == 3:
            self._sendSCPI(b"TRAC1:DEL:ALL \n")
            self._sendSCPI(b"TRAC2:DEL:ALL \n")
        
        elif (ch == 1) or (ch == 2):
            cmd = "TRAC{}:DEL:ALL \n".format(ch)
            self._sendSCPI(cmd.encode())
        
        else:
            raise ValueError("Channel must be one of 1,2,3(both).")

        pass

    def defineWaveform(
            self, id: int, length: int, init_value: int = 0,
            ch: int = 3) -> None:

        if length % 64 != 0:
            raise ValueError(
                "The length of a waveform must be a multiple of the "
                "gradularity.")

        if ch == 3:
            cmd = "TRAC{}:DEF {},{},{} \n".format(1, id, length, init_value)
            self._sendSCPI(cmd.encode())
            cmd = "TRAC{}:DEF {},{},{} \n".format(2, id, length, init_value)
            self._sendSCPI(cmd.encode())
        elif (ch == 1) or (ch == 2):
            cmd = "TRAC{}:DEF {},{},{} \n".format(ch, id, length, init_value)
            self._sendSCPI(cmd.encode())
        
    def setWaveform(
            self, string: str, id: int, ch: int = 3, offset: int = 0) -> None:

        if ch == 3:
            cmd = "TRAC{}:DATA {},{},{} \n".format(1, id, offset, string)
            self._sendSCPI(cmd.encode())
            cmd = "TRAC{}:DATA {},{},{} \n".format(2, id, offset, string)
            self._sendSCPI(cmd.encode())
        elif (ch == 1) or (ch == 2):
            cmd = "TRAC{}:DATA {},{},{} \n".format(ch, id, offset, string)
            self._sendSCPI(cmd.encode())

    def setFunctionMode(self, mode: str, ch: int = 3):
        options = {"Arbitary": "ARB", "Sequence": "STS", "Scenario": "STSC"}

        if mode not in options:
            raise ValueError(f"Mode must be one of: {options}. Not:{mode}.")
        
        if ch == 3:
            cmd = "FUNC1:MODE {} \n".format(options[mode])
            self._sendSCPI(cmd.encode())
            cmd = "FUNC2:MODE {} \n".format(options[mode])
            self._sendSCPI(cmd.encode())
        
        elif (ch == 1) or (ch == 2):
            cmd = "FUNC{}:MODE {} \n".format(ch, options[mode])
            self._sendSCPI(cmd.encode())
        
        else:
            raise ValueError("Channel must be one of 1,2,3(both).")

        pass

    def getFunctionMode(self, ch: int = 3):
        # options = {"ARB": "Arbitary", "STS": "Sequence", "STSC": "Scenario"}
        
        if ch == 3:
            mode1 = self._sendrecvSCPI(b"FUNC1:MODE? \n", 512)
            mode2 = self._sendrecvSCPI(b"FUNC2:MODE? \n", 512)

            return mode1, mode2
        
        elif (ch == 1) or (ch == 2):
            cmd = "FUNC{}:MODE? \n".format(ch)
            mode = self._sendrecvSCPI(cmd.encode(), 512)
            return mode
        
        else:
            raise ValueError("Channel must be one of 1,2,3(both).")

    def resetSequenceTable(self, ch: int = 3):
        
        if ch == 3:
            self._sendSCPI(b":STAB1:RES \n")
            self._sendSCPI(b":STAB2:RES \n")

        elif (ch == 1) or (ch == 2):
            cmd = ":STAB{}:RES \n".format(ch)
            self._sendSCPI(cmd.encode())

    def getSequenceTable(self, length: int, ch: int = 3):

        print(
            "WARNING: Due to AWG API limitations this should be read with a "
            "grain of salt.")

        length = int(length * 6)
        num_bytes = length * 11
        if ch == 3:
            seq_table1 = SequenceTable(self)
            seq_table2 = SequenceTable(self)
            cmd = "STAB1:DATA? 0,{} \n".format(length)
            cat1 = self._sendrecvSCPI(cmd.encode(), num_bytes)
            cmd = "STAB2:DATA? 0,{} \n".format(length)
            cat2 = self._sendrecvSCPI(cmd.encode(), num_bytes)
            seq_table1.read_from_str(cat1)
            seq_table2.read_from_str(cat2)

            return [seq_table1, seq_table2]

        elif (ch == 1) or (ch == 2):
            seq_table = SequenceTable(self)

            cmd = "STAB{}:DATA? 0,{} \n".format(int(ch), length)
            cat1 = self._sendrecvSCPI(cmd.encode(), num_bytes)
            seq_table.read_from_str(cat1)

            return seq_table

        else:
            raise ValueError("Channel must be one of 1,2,3(both).")

    def setSequenceTable(self, seqtable: SequenceTable, id: int, ch: int = 3):

        string = seqtable.write_to_str()

        cmd = "STAB{}:DATA {},{} \n"
        if ch == 3:
            self._sendSCPI((cmd.format(1, int(id), string)).encode())
            self._sendSCPI((cmd.format(2, int(id), string)).encode())

        elif (ch == 1) or (ch == 2):
            self._sendSCPI((cmd.format(ch, int(id), string)).encode())

    
class WaveformError(Exception):
    pass


class AWGError(Exception):
    pass


class SequenceError(Exception):
    pass
