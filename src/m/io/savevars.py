import shelve
import os.path


def savevars(*args):
    """
    SAVEVARS - function to save variables to a file.

    This function saves one or more variables to a file.  The names of the variables
    must be supplied.  If more than one variable is specified, it may be done with
    lists of names and values or a dictionary of name:value pairs.  All the variables
    in the workspace may be saved by specifying the globals() dictionary, but this
    may include a lot of extraneous data.

    Usage:
       savevars('shelve.dat', 'a', a)
       savevars('shelve.dat', ['a', 'b'], [a, b])
       savevars('shelve.dat', {'a':a, 'b':b})
       savevars('shelve.dat', globals())

    """

    filename = ''
    nvdict = {}

    if len(args) >= 1 and isinstance(args[0], str):
        filename = args[0]
        if not filename:
            filename = '/tmp/shelve.dat'

    else:
        raise TypeError("Missing file name.")

    if len(args) >= 3 and isinstance(args[1], str):  # (filename, name, value)
        for i in range(1, len(args), 2):
            nvdict[args[i]] = args[i + 1]

    elif len(args) == 3 and isinstance(args[1], list) and isinstance(args[2], list):  # (filename, [names], [values])
        for name, value in zip(args[1], args[2]):
            nvdict[name] = value

    elif len(args) == 2 and isinstance(args[1], dict):  # (filename, {names:values})
        nvdict = args[1]

    else:
        raise TypeError("Unrecognized input arguments.")

    if os.path.exists(filename):
        print(("Shelving variables to existing file '%s'." % filename))
    else:
        print(("Shelving variables to new file '%s'." % filename))

    my_shelf = shelve.open(filename, 'c')  # 'c' for create if not exist, else 'n' for new

    for name, value in list(nvdict.items()):
        try:
            my_shelf[name] = value
            print(("Variable '%s' shelved." % name))
        except TypeError:
            print(("Variable '%s' not shelved." % name))

    my_shelf.close()
