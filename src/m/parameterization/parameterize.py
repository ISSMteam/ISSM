import os
import datetime


def parameterize(md, parametername):
    """parameterize - parameterize a model

    From a parameter Python file, start filling in all the model fields that 
    were not filled in by the mesh.py and mask.py model methods. 
    Warning: the parameter file must be able to be run in Python

    Usage:
        md = parameterize(md, parametername)

    Example:
        md = parameterize(md, 'Square.py')
    """

    # Some checks
    if not os.path.exists(parametername):
        raise IOError("parameterize error message: file '%s' not found!" % parametername)

    # Try and run parameter file
    exec(compile(open(parametername).read(), parametername, 'exec'))

    # Name and notes
    if not md.miscellaneous.name:
        md.miscellaneous.name = os.path.basename(parametername).split('.')[0]

    md.miscellaneous.notes = 'Model created by using parameter file: \'%s\' on: %s.' % (parametername, datetime.datetime.strftime(datetime.datetime.now(), '%c'))

    return md
