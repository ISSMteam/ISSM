import numpy as np


def writejsfield(fid, name, variable, nods):
    #WRITEJSFIELD - write variable to javascript file
    #
    #   Usage:
    #      writejsfield(fid, name, variable)
    #
    #write array:
    #if not isinstance(variable, list):
    if type(variable[0]) == np.float64:
        fid.write('<!-- {0}{{{{{{-->\n'.format(name))
        fid.write('{0}=['.format(name))
        for i in range(0, nods - 1):
            fid.write('{0}, '.format(variable[i]))
        fid.write('{0}];\n'.format(variable[-1]))
        fid.write('<!--}}}}}}-->\n')
    else:
        #multi - sized array:
        fid.write('<!-- {0}{{{{{{-->\n'.format(name))
        fid.write('{0} = []\n'.format(name))
        for i in range(0, len(variable[2])):
            fid.write('{0}["{1}"] = ['.format(name, i))
            for j in range(1, nods - 1):
                fid.write('{0}, '.format(variable[j][i]))
            fid.write('{0}];\n'.format(variable[-1][i]))
        fid.write('<!--}}}}}}-->\n')
