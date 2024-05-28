from copy import deepcopy


class results(object):  #{{{
    """RESULTS class definition

    Usage:
        md.results = results()

    TODO:
    - Rework so that a solution of arbitrary length (normal or transient) can
    initialized from one call to the results class constructor.
    """

    def __init__(self):  #{{{
        pass
    # }}}

    def __repr__(self):  #{{{
        s = ''
        for key, value in self.__dict__.items():
            if isinstance(value, resultsdakota):
                lengthvalue = 1
            else:
                try:
                    lengthvalue = len(value)
                except TypeError:
                    lengthvalue = 1
            s += '    {}: [1x{} struct]\n'.format(key, lengthvalue)

        return s
    # }}}

    def setdefaultparameters(self):  #{{{
        #do nothing
        return self
    # }}}

    def checkconsistency(self, md, solution, analyses):  #{{{
        return md
    # }}}

    def marshall(self, prefix, md, fid):  #{{{
        pass
    # }}}
# }}}


class resultsdakota(object):  #{{{
    """RESULTSDAKOTA class definition - Used to store results from a run of
    Dakota.

    Usage:
        md.results.dakota = resultsdakota()

    NOTE: Values of attributes can themselves be instances of solution class.
    """

    def __init__(self):  #{{{
        pass
    # }}}

    def __repr__(self):  #{{{
        s = ''
        for key, value in self.__dict__.items():
            s += '    {}: '.format(key)
            if isinstance(value, list):
                s += '[{} element list]'.format(len(value))
            else:
                s += '{}'.format(value)
            s += '\n'
        return s
    # }}}

    def __len__(self):  #{{{
        return len(self.__dict__.keys())
    # }}}

    def setdefaultparameters(self):  #{{{
        #do nothing
        return self
    # }}}

    def checkconsistency(self, md, solution, analyses):  #{{{
        return md
    # }}}

    def marshall(self, prefix, md, fid):  #{{{
        pass
    # }}}
# }}}


class solution(object):  #{{{
    """SOLUTION class definition - Value of an attribute (which should be
    a string representing a solution type) of an instance of results class

    Elements of self.steps should be instances of solutionstep.

    Usage:
        setattr(md.results, 'SolutionType', solution())

    NOTE:
    - Under MATLAB, this is implemented as a structure array
    - Only when this instance of solution represents a transient solution
    should self.steps have more than one element
    """

    def __init__(self, *args):  #{{{
        self.steps = None
        if len(args) == 1:
            arg = args[0]
            if isinstance(arg, list):
                self.steps = arg
            else:
                raise Exception('solution class error: if initializing with an argument, that argument should be an empty list or a list of instances of solutionstep')
        else:
            self.steps = [solutionstep()]
    # }}}

    def __deepcopy__(self, memo):  #{{{
        return solution(deepcopy(self.steps, memo))
    # }}}

    def __repr__(self):  #{{{
        s = ''
        numsteps = len(self.steps)
        if numsteps == 1:
            for key, value in self.steps[0].__dict__.items():
                s += '    {}: {}\n'.format(key, value)
        else:
            s = '  1x{} struct array with fields:\n'.format(numsteps)
            s += '\n'
            for fieldname in self.steps[0].getfieldnames():
                s += '    {}\n'.format(fieldname)

        return s
    # }}}

    def __len__(self):  #{{{
        return len(self.steps)
    # }}}

    def __getattr__(self, key):  #{{{
        # NOTE: Currently only returning value from first frame of transient solution (see retrieval of md.results.TransientSolution.BedGRD in test2051.py for justification)
        return getattr(self.steps[0], key)

        # Original code
        # if len(self.steps) == 1:
        #     return getattr(self.steps[0], key)
        # else:
        #     raise Exception('<results>.<solution> error: Currently, can only get a field if we are not working with a transient solution.')
    # }}}

    def __getitem__(self, index):  #{{{
        while True:
            try:
                return self.steps[index]
            except:
                self.steps.append(solutionstep())
        else:
            raise Exception('<results>.<solution>: either request a specific result by index or make sure that there is only a single result for this solution (cannot be a transient solution)')
    # }}}

    def setdefaultparameters(self):  #{{{
        return self
    # }}}

    def checkconsistency(self, md, solution, analyses):  #{{{
        return md
    # }}}

    def marshall(self, prefix, md, fid):  #{{{
        pass
    # }}}
# }}}


class solutionstep(object):  #{{{
    """SOLUTIONSTEP class definition - Single element of <solution>.steps

    Usage:
        <solution>.steps.append(solutionstep())
    """

    def __init__(self, *args):  #{{{
        pass
    # }}}

    def __repr__(self):  #{{{
        s = ''
        width = self.getlongestfieldname()
        for key, value in self.__dict__.items():
            s += '    {:{width}s}: {}\n'.format(key, value, width=width)

        return s
    # }}}

    def getfieldnames(self):  #{{{
        return self.__dict__.keys()
    # }}}

    def getlongestfieldname(self):  #{{{
        maxlength = 0
        for key in self.__dict__.keys():
            length = len(key)
            if length > maxlength:
                maxlength = length

        return maxlength
    # }}}

    def setdefaultparameters(self):  #{{{
        return self
    # }}}

    def checkconsistency(self, md, solution, analyses):  #{{{
        return md
    # }}}

    def marshall(self, prefix, md, fid):  #{{{
        pass
    # }}}
# }}}
