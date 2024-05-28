from collections import OrderedDict


class pairoptions(object):
    """PAIROPTIONS class definition

    Usage:
        pairoptions = pairoptions()
        pairoptions = pairoptions('module', true, 'solver', false)
    """

    def __init__(self, *arg):  # {{{
        #self.functionname = ''
        self.list = OrderedDict()

        #get calling function name
        #import inspect
        #if len(inspect.stack()) > 1:
        #self.functionname = inspect.stack()[1][3]
        import traceback
        try:
            self.functionname = traceback.extract_stack(limit=2)[0][2]
        except IndexError:
            pass  #this is probably similar to the previous if treatment (but faster)

        #initialize list
        if not len(arg):
            pass  #Do nothing,
        else:
            self.buildlist(*arg)
    # }}}

    def __repr__(self):  # {{{
        s = "   functionname: '{}'\n".format(self.functionname)
        if self.list:
            s += "   list: ({}x{}) \n\n".format(len(self.list), 2)
            for item in self.list.items():
                s += "     field: {} value: '{}'\n".format(item[0], item[1])
            print(s)
        else:
            s += "   list: empty\n"
        return s
    # }}}

    def buildlist(self, *arg):  # {{{
        """BUILDLIST - build list of objects from input
        """

        #check length of input
        if len(arg) % 2:
            raise TypeError('Invalid parameter/value pair arguments')
        numoptions = int(len(arg) / 2)

        #go through arg and build list of objects
        for i in range(numoptions):
            if isinstance(arg[2 * i], str):
                self.list[arg[2 * i]] = arg[2 * i + 1]
            else:
                #option is not a string, ignore it
                print(("WARNING: option number {} is not a string and will be ignored.".format(i + 1)))
    # }}}

    def addfield(self, field, value):  # {{{
        """ADDFIELD - add a field to an options list
        """

        if isinstance(field, str):
            if field in self.list:
                print(("WARNING: field '{}' with value={} exists and will be overwritten with value={}.".format(field, str(self.list[field]), str(value))))
            self.list[field] = value
    # }}}

    def addfielddefault(self, field, value):  # {{{
        """ADDFIELDDEFAULT - add a field to an options list if it does not already exist
        """

        if isinstance(field, str):
            if field not in self.list:
                self.list[field] = value
    # }}}

    def AssignObjectFields(self, obj2):  # {{{
        """ASSIGNOBJECTFIELDS - assign object fields from options
        """

        for item in list(self.list.items()):
            if item[0] in dir(obj2):
                setattr(obj2, item[0], item[1])
            else:
                print(("WARNING: field '%s' is not a property of '%s'." % (item[0], type(obj2))))
        return obj2
    # }}}

    def changefieldvalue(self, field, newvalue):  # {{{
        """CHANGEOPTIONVALUE - change the value of an option in an option list
        """

        self.list[field] = newvalue
    # }}}

    def displayunused(self):  # {{{
        """DISPLAYUNUSED - display unused options
        """

        print('WARNING: pairoptions::displayunused is not yet implemented')
    # }}}

    def exist(self, field):  # {{{
        """EXIST - check if the option exists
        """

        #some argument checking:
        if field is None or field == '':
            raise ValueError('exist error message: bad usage')
        if not isinstance(field, str):
            raise TypeError("exist error message: field '%s' should be a string." % str(field))

        #Recover option
        if field in self.list:
            return True
        else:
            return False
    # }}}

    def getfieldvalue(self, field, default=None):  # {{{
        """GETFIELDVALUE - get the value of an option

        Usage:
           value = options.getfieldvalue(field, default)

        Find an option value from a field. A default option
        can be given in input if the field does not exist

        Examples:
           value = options.getfieldvalue(options, 'caxis')
           value = options.getfieldvalue(options, 'caxis', [0 2])
        """

        #some argument checking:
        if field is None or field == '':
            raise ValueError('getfieldvalue error message: bad usage')
        if not isinstance(field, str):
            raise TypeError("getfieldvalue error message: field '%s' should be a string." % str(field))

        #Recover option
        if field in self.list:
            value = self.list[field]
        else:
            if default is not None:
                value = default
            else:
                raise KeyError("error message: field '%s' has not been provided by user (and no default value has been specified)." % field)

        return value
    # }}}

    def removefield(self, field, warn):  # {{{
        """REMOVEFIELD - delete a field in an option list

        Usage:
           obj = removefield(self, field, warn)

        if warn == 1 display an info message to warn user that
        some of their options have been removed.
        """

        #check if field exist
        if field in self.list:

            #remove duplicates from the options list
            del self.list[field]

            #warn user if requested
            if warn:
                print(("removefield info: option '%s' has been removed from the list of options." % field))
    # }}}

    def marshall(self, md, fid, firstindex):  # {{{

        for i, item in enumerate(self.list.items()):
            name = item[0]
            value = item[1]

            raise NameError('need to sync with MATLAB')

        #Write option name
        #WriteData(fid, prefix, 'enum', (firstindex - 1) + 2 * i + 1, 'data', name, 'format', 'String')

        #Write option value
        #if   isinstance(value, (str, unicode)):
        #    WriteData(fid, prefix, 'enum', (firstindex - 1) + 2 * i + 2, 'data', value, 'format', 'String')
        #elif isinstance(value, (bool, int, long, float)):
        #    WriteData(fid, prefix, 'enum', (firstindex - 1) + 2 * i + 2, 'data', value, 'format', 'Double')
        #else:
        #raise TypeError("Cannot marshall option '%s': format not supported yet." % name)
    # }}}
