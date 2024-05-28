import numpy as np

from MatlabArray import string_dim


class histogram_bin_uncertain(object):
    '''
    HISTOGRAM_BIN_UNCERTAIN class definition

    Usage:
        [hbu] = histogram_bin_uncertain(
            'descriptor', descriptor,
            'pairs_per_variable', pairs_per_variable,
            'abscissas', abscissas,
            'counts', counts
            )

        where the required args are:
            descriptor          (char, description, '')
            pairs_per_variable  (double list, [])
            abscissas           (double list, [])
            counts              (int list, [])

    NOTE: A call to the constructor with zero arguments will return a default 
    instance; one argument of the class copies the instance; three or more 
    arguments constructs a new instance from the arguments.
    '''

    def __init__(self):  # {{{
        self.descriptor = ''
        self.pairs_per_variable = []
        self.abscissas = []
        self.counts = []
    # }}}

    @staticmethod
    def histogram_bin_uncertain(*args):  # {{{
        nargin = len(args)

        # create a default object
        if nargin == 0:
            return histogram_bin_uncertain()

        # copy the object
        elif nargin == 1:
            if isinstance(args[0], histogram_bin_uncertain):
                hbu = args[0]
            else:
                raise Exception("Object {} is a {} class object, not 'histogram_bin_uncertain'.".format(str(args[0]), str(type(args[0]))))

        elif nargin == 2 or nargin == 3:
            raise Exception("Construction of 'histogram_bin_uncertain' class object requires at least {} inputs.".format(4))

        # create the object from the input
        elif nargin == 4:
            hbu = histogram_bin_uncertain()

            #recover options:
            options = pairoptions(*args)

            #initialize fields:
            hbu.descriptor          = options.getfieldvalue('descriptor')
            hbu.pairs_per_variable  = options.getfieldvalue('pairs_per_variable')
            hbu.abscissas           = options.getfieldvalue('abscissas')
            hbu.counts              = options.getfieldvalue('counts')

        else:
            raise Exception("Construction of histogram_bin_uncertain class object requires either (1) no arguments, (2) a histogram_bin_uncertain instance to copy from, or (3) a descriptor and pairs per variable, abscissas, and counts lists")

    @staticmethod
    def __repr__(hbu):  # {{{
        s = ""
        for i in range(len(hbu)):
            s += "class {} object {} = \n".format(hbu.__class__.__name__, string_dim(hbu, i))
        s = "{}\n{}".format(s, fielddisplay(self, 'descriptor', 'name tag'))
        s = "{}\n{}".format(s, fielddisplay(self, 'pairs_per_variable', 'pairs per variable'))
        s = "{}\n{}".format(s, fielddisplay(self, 'abscissas', 'abscissas'))
        s = "{}\n{}".format(s, fielddisplay(self, 'counts', 'counts'))

        return s
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        return
    # }}}

    #virtual functions needed by qmu processing algorithms
    #implemented:

    @staticmethod
    def prop_desc(hbu, dstr):  # {{{
        desc = ['' for i in range(np.size(hbu))]
        for i in range(np.size(hbu)):
            if hbu[i].descriptor != '' or type(hbu[i].descriptor) != str:
                desc[i] = str(hbu[i].descriptor)
            elif dstr != '':
                desc[i] = str(dstr) + str(string_dim(hbu, i, 'vector'))
            else:
                desc[i] = 'hbu' + str(string_dim(hbu, i, 'vector'))

        desc = allempty(desc)

        return desc
    # }}}

    @staticmethod
    def prop_mean(hbu):  # {{{
        mean = np.zeros(np.size(hbu))
        for i in range(np.size(hbu)):
            mean[i] = hbu[i].mean
        return mean
    # }}}

    @staticmethod
    def prop_stddev(hbu):  # {{{
        stddev = np.zeros(np.size(hbu))
        for i in range(np.size(hbu)):
            stddev[i] = hbu[i].stddev
        return stddev
    # }}}

    @staticmethod
    def prop_lower(hbu):  # {{{
        lower = []
        return
    # }}}

    @staticmethod
    def prop_upper(hbu):  # {{{
        upper = []
        return upper
    # }}}

    #default
    @staticmethod
    def prop_abscissas(hbu):  # {{{
        abscissas = []
        for i in range(len(hbu)):
            abscissas.extend(hbu[i].abscissas)
        abscissas = allequal(abscissas, -np.inf)
        return abscissas
    # }}}

    @staticmethod
    def prop_pairs_per_variable(hbu):  # {{{
        pairs_per_variable = np.zeros((1, len(hbu)))
        for i in range(len(hbu)):
            pairs_per_variable[i] = hbu[i].pairs_per_variable
        abscissas = allequal(pairs_per_variable, -np.inf)
        return pairs_per_variable
    # }}}

    @staticmethod
    def prop_counts(hbu):  # {{{
        counts = []
        for i in range(len(hbu)):
            counts.extend(hbu[i].counts)
        counts = allequal(counts, -np.inf)
        return counts
    # }}}

    @staticmethod
    def prop_initpt(hbu):  # {{{
        initpt = []
        return initpt
    # }}}

    @staticmethod
    def prop_initst(hbu):  # {{{
        inist = []
        return inist
    # }}}

    @staticmethod
    def prop_stype(hbu):  # {{{
        stype = []
        return stype
    # }}}

    @staticmethod
    def prop_scale(hbu):  # {{{
        scale = []
        return scale
    # }}}

    #new methods:
    def isscaled(self):  # {{{
        if strncmp(self.descriptor, 'scaled_', 7):
            return True
        else:
            return False
    # }}}

    @staticmethod
    def dakota_write(fidi, dvar):
        # possible namespace pollution, the above import seems not to work
        from vlist_write import vlist_write
        # collect only the variables of the appropriate class
        hbu = deepcopy(dvar)
        fields = fieldnames(hbu)
        for field in fields:
            if getattr(hbu, field)[0].__class__.__name__ != 'histogram_bin_uncertain':
                delattr(hbu, field)
        if len(hbu) > 0:
            vlist_write(fidi, 'histogram_bin_uncertain', 'hbu', hbu)
    # }}}
