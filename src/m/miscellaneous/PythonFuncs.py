import numpy as np


def logical_and_n(*arg):  # {{{
    if len(arg):
        result = arg[0]
        for item in arg[1:]:
            result = np.logical_and(result, item)
        return result
    else:
        return None
# }}}

def logical_or_n(*arg):  # {{{
    if len(arg):
        result = arg[0]
        for item in arg[1:]:
            result = np.logical_or(result, item)
        return result
    else:
        return None
# }}}
