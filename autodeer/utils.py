import re

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



