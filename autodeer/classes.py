import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import hilbert
from scipy.io import savemat
import os
from autodeer import __version__
import copy
import time
import datetime
import numbers
import uuid
import json
import base64
from autodeer.utils import autoEPRDecoder
from pathlib import Path
import logging

# =============================================================================


class Interface:
    """Represents the interface connection from autoEPR to the spectrometer.
    """

    def __init__(self,log=None) -> None:
        self.pulses = {}
        self._savefolder = str(Path.home())
        self.savename = ""
        if log is None:
            self.log = logging.getLogger('interface')
        else:
            self.log = log
        pass

    def connect(self) -> None:
        pass

    def acquire_dataset(self, data):
        """
        Acquires the dataset.
        """

        # data.sequence = self.cur_exp
        data.attrs['time'] = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        return data
    

    def launch(self, sequence, savename: str):
        """Launches the experiment and initialises autosaving.

        Parameters
        ----------
        sequence : Sequence
            The sequence to be launched
        savename : str
            The savename for this measurement. A timestamp will be added to the value.
        """
        timestamp = datetime.datetime.now().strftime(r'%Y%m%d_%H%M_')
        self.savename=timestamp + savename + '.h5'
        pass

    def isrunning(self) -> bool:
        return False

    def terminate(self) -> None:
        """
        Terminates the experiment immediately. 
        """
        pass

    def terminate_at(self, criterion, test_interval=2, keep_running=True, verbosity=0,autosave=True):
        """Terminates the experiment upon a specific condition being
        satisified. 

        Parameters
        ----------
        criterion : _type_
            The criteria to be tested.
        test_interval : int, optional
            How often should the criteria be tested in minutes, by default 10.
        keep_running : bool, optional
            If True, an error will not be raised if the experiment finishes before the criteria is met, by default True.
        verbosity : int, optional
            The verbosity level, by default 0. 
        autosave : bool, optional
            If True, the data will be autosaved, by default True.

        """



        test_interval_seconds = test_interval * 60
        condition = False
        last_scan = 0


        while not condition:

            if not self.isrunning():
                if keep_running:
                    self.terminate()
                    return None
                else:
                    msg = "Experiments has finished before criteria met."
                    raise RuntimeError(msg)

            start_time = time.time()
            data = self.acquire_dataset()
            if autosave:
                self.log.debug(f"Autosaving to {os.path.join(self._savefolder,self.savename)}")
                data.to_netcdf(os.path.join(self._savefolder,self.savename),engine='h5netcdf',invalid_netcdf=True)

            try:
                # nAvgs = data.num_scans.value
                nAvgs = data.attrs['nAvgs']
            except AttributeError:
                self.log.warning("WARNING: Dataset missing number of averages(nAvgs)!")
                nAvgs = 1
            finally:
                if nAvgs < 1:
                    time.sleep(30)  # Replace with single scan time
                    continue
                elif nAvgs <= last_scan:
                    time.sleep(30)
                    continue    
            last_scan = nAvgs
            if verbosity > 0:
                print("Testing")
            condition = criterion.test(data, verbosity)

            if not condition:
                end_time = time.time()
                if (end_time - start_time) < test_interval_seconds:
                    if verbosity > 0:
                        print("Sleeping")
                    time.sleep(test_interval_seconds - (end_time - start_time))

        if callable(criterion.end_signal):
            criterion.end_signal()
        self.terminate()
        pass

# =============================================================================

class Parameter:
    """
    Represents a sequence or pulse parameter.
    """

    def __init__(self, name, value, unit="", description="", virtual=False,
                  **kwargs) -> None:
        """A general parameter.

        Parameters
        ----------
        name : str
            The parameter name
        value : float or int
            The parameter value, eithe initial or static
        unit : str, optional
            The unit of parameter, by default None. Leave as None if unitless.
        description : str, optional
            A brief description of the parameter, by default None    
        axis : np.ndarray, optional
            The difference from the intial value for each position in a 
            dynamic axis. Can be n-dimensional, by default None.
        ax_id : list, optional 
        virtual: bool, optional
            A virtual paramter is only used to vary other parameters, it is not 
            varied itself and will not be directly passed to a spectrometer.
            This parameter is **never** inherited.
            By default, False
             
    

        Attributes
        ----------
        progressive : bool
            Is the parameter used in any progression or is it constant
        prog : dict
            A dict containing progressive programs for this parameter. This 
            list has two elements. 1) The axis_id"s and 2) the "axis" of 
            values.

        Parameter Arthimatic
        --------------------



        Examples
        --------
        Creating a static parameter
        ```
        Par1 = Parameter(
            name="Par1", value=10, unit="us", description="The first parameter")
        ```
        Creating a dynamic parameter
        ```
        Par1 = Parameter(
            name="Par1", value=10, unit="us", description="The first parameter",
            axis=np.arange(0,10,1), axis_id=0)
        ```

        Adding a parameter and a number:
        ```
        Par1 = Parameter(
            name="Par1", value=10, unit="us", description="The first parameter")
        Par2 = Par1 + 2
        """

        self.name = name

        if isinstance(value, Parameter):
            self.value = value.value
        else:
            self.value = value
            
        self.unit = unit
        self.description = description
        self.virtual = virtual
        self.axis = []
        self.ax_id = []
        if "link" in kwargs:
            if not isinstance(kwargs["link"], Parameter):
                raise ValueError("The linked parameter must be a Parmater object")
            self.uuid = kwargs["link"].uuid
        else:
            self.uuid = uuid.uuid1()

        if "step" in kwargs:
            step = kwargs["step"]
            dim = kwargs["dim"]
            if "axis_id" in kwargs:
                axis_id = kwargs["axis_id"]
            else:
                axis_id = 0
            if "start" in kwargs:
                start = kwargs["start"]
            else:
                start = 0
            if step == 0:
                axis = np.zeros(dim)
            else:
                axis = np.arange(start=start, stop= dim*step+start,step=step)
            self.add_axis(axis=axis,axis_id=axis_id)
        


        pass

    def add_axis(self, axis_id, axis):
        # if self.axis == []:
        #     self.axis.append(np.array(axis))
        #     self.ax_id.append(axis_id)
        self.axis.append({"axis":axis, "uuid":self.uuid})

    def get_axis(self):
        init_value = self.value
        axes = []
        for axis in self.axis:
            axes.append(axis['axis'] + init_value)
        if len(axes) == 1:
            return axes[0]
        else:
            return axes
    
    @property
    def dim(self):
        if self.axis is []:
            return ()
        dims = []
        for ax in self.axis:
            dims.append(ax['axis'].shape[0])
        return tuple(dims)

    def remove_dynamic(self):
        self.axis = []
        self.ax_id =[]     
    
    def is_static(self) -> bool:
        if self.axis == []:
            return True
        # elif self.ax_id == []:
        #     return True
        else:
            return False

    def __eq__(self, __o: object) -> bool:
        if type(__o) is not Parameter:
            raise ValueError(
                "Equivalence only works between Parameter classes")
        return self.value == __o.value
    
    def __add__(self, __o:object):

        if type(__o) is Parameter:
            if self.unit != __o.unit:
                raise RuntimeError("Both parameters must have the same unit")
            new_value = self.value + __o.value
            new_name = f"{self.name} + {__o.name}"
            new_description  = new_name
            new_parameter = Parameter(
                name=new_name, value=new_value, unit=self.unit,
                description=new_description)
            if not self.is_static():
                if not __o.is_static():
                    # Dynamic parmaters can only be summed and multiplied if the axis has the same uuid. I.e. they were linked when created or are deriratives of each other. 
                    new_ax_id = []
                    new_axis = []
                    # a_ax_ids:list = self.ax_id
                    a_ax_ids:list = [self.axis[i]["uuid"] for i in range(len(self.axis))]
                    # b_ax_ids:list = __o.ax_id
                    b_ax_ids:list = [__o.axis[i]["uuid"] for i in range(len(__o.axis))]
                    ab_ax_ids = list(set(a_ax_ids + b_ax_ids))
                    for id in ab_ax_ids:
                        if id not in b_ax_ids: # I.e. only in A
                            a_index = a_ax_ids.index(id)
                            new_axis.append(self.axis[a_index])
                            new_ax_id.append(id)
                        elif id not in a_ax_ids: # I.e. only in B
                            b_index = b_ax_ids.index(id)
                            new_axis.append(__o.axis[b_index])
                            new_ax_id.append(id)
                        else: # in both
                            a_index = a_ax_ids.index(id)
                            b_index = b_ax_ids.index(id)
                            b_ax_ids.remove(id)
                            new_axis.append({"axis": self.axis[a_index]["axis"] + __o.axis[b_index]["axis"], "uuid": id})
                            new_ax_id.append(id)
                else:
                    new_axis = self.axis
                    new_ax_id = self.ax_id

            else:
                if not __o.is_static():
                    new_axis = __o.axis
                    new_ax_id = __o.ax_id
                else:
                    new_axis = []
                    new_ax_id = []

            new_parameter.axis = new_axis
            new_parameter.ax_id = new_ax_id

            return new_parameter

        elif isinstance(__o, numbers.Number):
            new_value = self.value + __o
            new_name = f"{self.name} + {__o}"
            new_parameter = Parameter(
                name=new_name, value=new_value, unit=self.unit)
            if not self.is_static():
                new_axis = self.axis
                new_ax_id = self.ax_id
                new_parameter.axis = new_axis
                new_parameter.ax_id = new_ax_id
            return new_parameter
        
        elif isinstance(__o, np.ndarray):
            if self.axis.shape != __o.shape:
                raise RuntimeError(
                    "Both parameters axis and the array must have the same shape")
    
    def __sub__(self, __o:object):
        
        if type(__o) is Parameter:
            if self.unit != __o.unit:
                raise RuntimeError("Both parameters must have the same unit")
            new_value = self.value - __o.value
            new_name = f"{self.name} - {__o.name}"
            new_parameter = Parameter(
                name=new_name, value=new_value, unit=self.unit)
            if not self.is_static():
                if not __o.is_static():
                    # Dynamic parmaters can only be summed and multiplied if the axis has the same uuid. I.e. they were linked when created or are deriratives of each other. 
                    new_ax_id = []
                    new_axis = []
                    # a_ax_ids:list = self.ax_id
                    a_ax_ids:list = [self.axis[i]["uuid"] for i in range(len(self.axis))]
                    # b_ax_ids:list = __o.ax_id
                    b_ax_ids:list = [__o.axis[i]["uuid"] for i in range(len(self.__o))]
                    ab_ax_ids = list(set(a_ax_ids + b_ax_ids))
                    for id in ab_ax_ids:
                        if id not in b_ax_ids: # I.e. only in A
                            a_index = a_ax_ids.index(id)
                            new_axis.append({"axis": self.axis[a_index], "uuid": self.uuid})
                            new_ax_id.append(id)
                        elif id not in a_ax_ids: # I.e. only in B
                            b_index = b_ax_ids.index(id)
                            new_axis.append({"axis": __o.axis[b_index], "uuid": __o.uuid})
                            new_ax_id.append(id)
                        else: # in both
                            a_index = a_ax_ids.index(id)
                            b_index = b_ax_ids.index(id)
                            b_ax_ids.remove(id)
                            new_axis.append({"axis": self.axis[a_index] - __o.axis[b_index], "uuid": id})
                            new_ax_id.append(id)
                else:
                    new_axis = self.axis
                    new_ax_id = self.ax_id

            else:
                if not __o.is_static():
                    new_axis = __o.axis
                    new_ax_id = __o.ax_id
                else:
                    new_axis = []
                    new_ax_id = []

            new_parameter.axis = new_axis
            new_parameter.ax_id = new_ax_id

            return new_parameter

        elif isinstance(__o, numbers.Number):
            new_value = self.value - __o
            new_name = f"{self.name} - {__o}"
            new_parameter = Parameter(
                name=new_name, value=new_value, unit=self.unit)
            if self.axis is not None:
                new_axis = self.axis 
                new_ax_id = self.ax_id
                new_parameter.axis = new_axis
                new_parameter.ax_id = new_ax_id
            return new_parameter
        
        elif isinstance(__o, np.ndarray):
            if self.axis.shape != __o.shape:
                raise RuntimeError(
                    "Both parameters axis and the array must have the same shape")

    def __mul__(self, __o:object):
        if type(__o) is Parameter:
            if self.unit != __o.unit:
                raise RuntimeError("Both parameters must have the same unit")
            # if not __o.is_static():
            #     raise RuntimeError("Multiplictaion of two dynamic parameters is not supported")
            new_value = self.value * __o.value
            new_name = f"{self.name} * {__o.name}"
            new_parameter = Parameter(
                name=new_name, value=new_value, unit=self.unit)
            # if self.axis is not None:
            #     new_axis =  [np.array([item * __o.value for item in axis]) for axis in self.axis ]
            #     new_ax_id = self.ax_id
            #     new_parameter.axis = new_axis
            #     new_parameter.ax_id = new_ax_id
            # return new_parameter
            if not self.is_static():
                if not __o.is_static():
                    # Dynamic parmaters can only be summed and multiplied if the axis has the same uuid. I.e. they were linked when created or are deriratives of each other. 
                    new_ax_id = []
                    new_axis = []
                    # a_ax_ids:list = self.ax_id
                    a_ax_ids:list = [self.axis[i]["uuid"] for i in range(len(self.axis))]
                    # b_ax_ids:list = __o.ax_id
                    b_ax_ids:list = [__o.axis[i]["uuid"] for i in range(len(self.__o))]
                    ab_ax_ids = list(set(a_ax_ids + b_ax_ids))
                    for id in ab_ax_ids:
                        if id not in b_ax_ids: # I.e. only in A
                            a_index = a_ax_ids.index(id)
                            new_axis.append({"axis": self.axis[a_index], "uuid": self.uuid})
                            new_ax_id.append(id)
                        elif id not in a_ax_ids: # I.e. only in B
                            b_index = b_ax_ids.index(id)
                            new_axis.append({"axis": __o.axis[b_index], "uuid": __o.uuid})
                            new_ax_id.append(id)
                        else: # in both
                            a_index = a_ax_ids.index(id)
                            b_index = b_ax_ids.index(id)
                            b_ax_ids.remove(id)
                            new_axis.append({"axis": self.axis[a_index] * __o.axis[b_index], "uuid": id})
                            new_ax_id.append(id)
                else:
                    new_axis = self.axis
                    new_ax_id = self.ax_id

            else:
                if not __o.is_static():
                    new_axis = __o.axis
                    new_ax_id = __o.ax_id
                else:
                    new_axis = []
                    new_ax_id = []

            new_parameter.axis = new_axis
            new_parameter.ax_id = new_ax_id

            return new_parameter

        elif isinstance(__o, numbers.Number):
            new_value = self.value * __o
            new_name = f"{self.name} + {__o}"
            new_parameter = Parameter(
                name=new_name, value=new_value, unit=self.unit)
            if self.axis is not []:
                new_axis = copy.deepcopy(self.axis)
                for i,axis in enumerate(new_axis):
                    new_axis[i]["axis"] = np.array(axis["axis"])* __o
                # new_axis = [np.array([item * __o for item in axis]) for axis in self.axis ]
                new_ax_id = self.ax_id
                new_parameter.axis = new_axis
                new_parameter.ax_id = new_ax_id
            return new_parameter
        
        elif isinstance(__o, np.ndarray):
            if self.axis.shape != __o.shape:
                raise RuntimeError(
                    "Both parameters axis and the array must have the same shape")

    def __rmul__(self, __o:object):
        return self.__mul__(__o)

    def copy(self):
        return copy.deepcopy(self)
            
    def _to_dict(self):
        to_return = {"version": __version__, "type": "Parameter"}


        for key, var in vars(self).items():
            if isinstance(var, np.ndarray):
                data_b64 = base64.b64encode(var.data)
                to_return[key] =  dict(__ndarray__=str(data_b64),
                            dtype=str(var.dtype),
                            shape=var.shape)
            if isinstance(var, complex):
                to_return[key] = str(var)
            if isinstance(var, uuid.UUID):
                return_dict = {"__uuid__": str(var)}
                to_return[key] = return_dict
            elif isinstance(var, Parameter):
                to_return[key] = var._to_dict()
            else:
                to_return[key] = var

        return to_return
    
    def _to_json(self):
        class autoEPREncoder(json.JSONEncoder):
            def default(self, obj):
                if isinstance(obj, np.ndarray):
                    if (len(obj) > 0 ) and isinstance(obj[0], str):
                        return list(obj)
                    data = np.ascontiguousarray(obj.data)
                    data_b64 = base64.b64encode(data)
                    return dict(__ndarray__=str(data_b64),
                                dtype=str(obj.dtype),
                                shape=obj.shape)
                if isinstance(obj, complex):
                    return str(obj)
                if isinstance(obj, numbers.Number):
                    return float(obj)
                if isinstance(obj, uuid.UUID):
                    return_dict = {"__uuid__": str(obj)}
                    return return_dict
                if isinstance(obj, Parameter):
                    return obj._to_dict()
                else:
                    return json.JSONEncoder.default(self, obj)
        
        return json.dumps(self._to_dict(), cls=autoEPREncoder, indent=4)
    
    def save(self, filename):
        """Save the parameter to a JSON file.

        Parameters
        ----------
        filename : str
            Path to the JSON file.

        Returns
        -------
        None

        Raises
        ------
        TypeError
            If the object cannot be serialized to JSON.

        Example
        -------
        >>> obj = Parameter()
        >>> obj.save("my_parameter.json")
        """
        with open(filename, "w") as f:
           f.write(self._to_json())

    @staticmethod
    def _from_dict(dict):
          new_param = Parameter(
              name=dict['name'], value=dict['value'], 
              unit=dict['unit'], description=dict['description'])
          new_param.axis = dict['axis']
          new_param.ax_id = dict['ax_id']
          new_param.uuid = dict['uuid']
          return new_param
    

    @classmethod
    def _from_json(cls, JSONstring):
        dct = json.loads(JSONstring, object_hook=autoEPRDecoder)
        return cls._from_dict(dct)
    
    @classmethod
    def load(cls, filename):
        """Load a Parameter object from a JSON file.

        Parameters
        ----------
        filename : str
            Path to the JSON file.

        Returns
        -------
        obj : Parameter
            The Pulse loaded from the JSON file.

        Raises
        ------
        FileNotFoundError
            If the file does not exist.

        Example
        -------
        >>> obj = Parameter.load("my_parameter.json")
        """
        with open(filename, "r") as f:
           file_buffer = f.read()
        return cls._from_json(file_buffer)
          
