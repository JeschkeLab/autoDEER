import re
import numpy as np
from scipy.sparse import bsr_array
import uuid
import base64


def build_table(source, params, params_widths):
    string = ""
    params_used = []
    title_str = ""
    line_str = ""
    title_fmt = []
    for i, param in enumerate(params):
        if hasattr(source[0], param):
            attr = getattr(source[0], param)
            line_str += f" {{:<{params_widths[i]}}}"
            tmp = re.findall(r'(\d+)', params_widths[i])[0]
            title_str += f" {{:<{tmp}}}"
            if attr.unit is None:
                title_fmt.append(attr.name)
            else:
                title_fmt.append(attr.name + ' (' + attr.unit + ")")
            params_used.append(param)
        elif param in ['iD', 'type' ,'Phase Cycle']:
            line_str += f" {{:<{params_widths[i]}}}"
            tmp = re.findall(r'(\d+)', params_widths[i])[0]
            title_str += f" {{:<{tmp}}}"
            title_fmt.append(param)
            params_used.append(param)

    title_str += "\n"
    line_str += "\n"
    string += title_str.format(*title_fmt)

    for i, pulse in enumerate(source):
        elements = []
        for param in params_used:
            if param == "type":
                elements.append(type(pulse).__name__)
            elif param == "iD":
                elements.append(i)
            elif param == 'Phase Cycle':
                elements.append(pulse._pcyc_str())
            elif hasattr(pulse, param):
                if getattr(pulse, param) is None:
                    elements.append("N/A")
                elif getattr(pulse, param).value is None:
                    elements.append("None")
                else:
                    elements.append(f"{getattr(pulse, param).value:>5.5g}")
            else:
                elements.append("N/A")
        string += line_str.format(*elements)

    return string


def sop(spins, comps):
    """Spin Operator Matricies.

    This function is ported from EasySpin (https://easyspin.org/easyspin/documentation/sop.html) 

    References:
    +++++++++++
    [1] Stefan Stoll, Arthur Schweiger
    EasySpin, a comprehensive software package for spectral simulation and analysis in EPR
    J. Magn. Reson. 178(1), 42-55 (2006)
    
    [2] Stefan Stoll, R. David Britt
    General and efficient simulation of pulse EPR spectra
    Phys. Chem. Chem. Phys. 11, 6614-6625 (2009)

    Parameters
    ----------
    spins : list
        A list of each spin and its spin qunatum number
    comps : str
        The type of spin operator matrix to create. Options are: x,y,z,+,-,e
    """
    
    num_spins = len(spins)
    OP=np.array([1])

    for spin_num in range(num_spins):
        I = spins[spin_num]
        sop_type = comps[spin_num]
        n = int(I * 2 + 1)
        if sop_type == 'x':
            m = np.arange(1,n)
            r = np.hstack((m-1,m))
            c = np.hstack((m,m-1))
            dia = 0.5 * np.sqrt(m*m[::-1])
            val = np.hstack((dia,dia))
        elif sop_type == 'y':
            m = np.arange(1,n)
            r = np.hstack((m-1,m))
            c = np.hstack((m,m-1))
            dia = -0.5*1j * np.sqrt(m*m[::-1])
            val = np.hstack((dia,-dia))
        elif sop_type == 'z':
            m = np.arange(1, n+ 1)
            r = m - 1
            c = m - 1
            val = -m + I + 1
        elif sop_type == '+':
            m = np.arange(1,n)
            r = m -1
            c = m
            val = np.sqrt(m * m[::-1])
        elif  sop_type == '-':
            m = np.arange(1,n)
            r = m + 1
            c = m
            val = np.sqrt(m * m[::-1])
        elif  sop_type == 'e':
            m = np.arange(1,n+1)
            r = m-1
            c = m-1
            val = np.ones(n)
        else:
            raise ValueError(f"Incorect specification of comps: ",
                             f"{sop_type} is not a valid input")
        M_ = bsr_array((val, (r,c)), shape=(n,n)).toarray()
        OP = np.kron(OP, M_)

    return OP

def transpose_dict_of_list(d):
    """Turns a dictionary of lists into a list of dictionaries.
    """
    return [dict(zip(d, col)) for col in zip(*d.values())]

def transpose_list_of_dicts(d):
    """Turns a list of dictionaries into a dictionary of lists.
    """

    if len(d) == 0:
        return {}
    else:
        return {key: [i[key] for i in d] for key in d[0]}

def save_file(path, str):
    with open(path, "w") as file:
        file.write(str)


def autoEPRDecoder(dct):
    if isinstance(dct, dict) and '__uuid__' in dct:
        return uuid.UUID(dct["__uuid__"])
    if isinstance(dct, dict) and '__ndarray__' in dct:
        data = base64.b64decode(dct['__ndarray__'][2:-1])
        return np.frombuffer(data, dct['dtype']).reshape(dct['shape'])
    return dct


def gcd(values:list):
    """Generates the greatest common dividor on a  list of floats

    Parameters
    ----------
    values : list
        _description_
    """

    if len(values) == 1:
        return values[0]
    
    if len(values) == 2:
        a = values[0]
        b = values[1]
        while b:
            a, b = b, a % b
        return a
    
    if len(values) > 2:
        a = values[0]
        b = values[1]
        while b:
            a, b = b, a % b
        values[0] = a
        values.pop(1)
        return gcd(values)
    

def val_in_us(Param):
    """Returns the value or axis of a parameter in microseconds

    Parameters
    ----------
    Param : autodeer.Parameter
        The parameter to be converted

    Returns
    -------
    float or numpy.ndarray
    """
    if len(Param.axis) == 0:
        if Param.unit == "us":
            return Param.value
        elif Param.unit == "ns":
            return Param.value / 1e3
    elif len(Param.axis) == 1:
        if Param.unit == "us":
            return Param.value + Param.axis[0]['axis']
        elif Param.unit == "ns":
            return (Param.value + Param.axis[0]['axis']) / 1e3 
    else:
        raise ValueError("Parameter must have 0 or 1 axes")

def val_in_ns(Param):
    """Returns the value or axis of a parameter in nanoseconds

    Parameters
    ----------
    Param : autodeer.Parameter
        The parameter to be converted

    Returns
    -------
    float or numpy.ndarray
    """
     
    if len(Param.axis) == 0:
        if Param.unit == "us":
            return Param.value * 1e3
        elif Param.unit == "ns":
            return Param.value 
    elif len(Param.axis) == 1:
        if Param.unit == "us":
            return (Param.tau1.value + Param.axis[0]['axis']) * 1e3
        elif Param.unit == "ns":
            return (Param.value + Param.axis[0]['axis']) 
    else:
        raise ValueError("Parameter must have 0 or 1 axes")

