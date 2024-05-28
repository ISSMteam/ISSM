from pairoptions import pairoptions
from WriteData import WriteData


class verbose(object):
    """VERBOSE class definition

    Available verbosity levels:
        mprocessor  : model processing
        module      : modules
        solution    : solution sequence
        solver      : solver info (extensive)
        convergence : convergence criteria
        control     : control method
        qmu         : sensitivity analysis
        autodiff    : AD analysis
        smb         : SMB analysis

    Usage:
        verbose = verbose()
        verbose = verbose(3)
        verbose = verbose('001100')
        verbose = verbose('module', True, 'solver', False)

    NOTE:
    - Some parts of this file are Synchronized with 
    src/c/shared/Numerics/Verbosity.h. Do not modify these sections. See 
    src/c/shared/Numerics/README for more info.
    """

    def __init__(self, *args):  # {{{
        #BEGINFIELDS
        self.mprocessor = False
        self.module = False
        self.solution = False
        self.solver = False
        self.convergence = False
        self.control = False
        self.qmu = False
        self.autodiff = False
        self.smb = False
        #ENDFIELDS

        if not len(args):
            # Don't do anything
            self.solution = True
            self.qmu = True
            self.control = True
            pass

        elif len(args) == 1:
            binary = args[0]
            if isinstance(binary, str):
                if binary.lower() == 'all':
                    binary = pow(2, 11) - 1  # all ones
                    self.BinaryToVerbose(binary)
                    self.solver = False  # Do not use by default
                else:
                    binary = int(binary, 2)
                    self.BinaryToVerbose(binary)
            elif isinstance(binary, (int, float)):
                self.BinaryToVerbose(int(binary))

        else:
            # Use options to initialize object
            self = pairoptions(*args).AssignObjectFields(self)

            # Cast to logicals
            listproperties = vars(self)
            for fieldname, fieldvalue in list(listproperties.items()):
                if isinstance(fieldvalue, bool) or isinstance(fieldvalue, (int, float)):
                    setattr(self, fieldname, bool(fieldvalue))
                else:
                    raise TypeError("verbose supported field values are logicals only (True or False)")
    # }}}

    def __repr__(self):  # {{{

        #BEGINDISP
        s = "class '%s' = \n" % type(self)
        s += "   %15s : %s\n" % ('mprocessor', self.mprocessor)
        s += "   %15s : %s\n" % ('module', self.module)
        s += "   %15s : %s\n" % ('solution', self.solution)
        s += "   %15s : %s\n" % ('solver', self.solver)
        s += "   %15s : %s\n" % ('convergence', self.convergence)
        s += "   %15s : %s\n" % ('control', self.control)
        s += "   %15s : %s\n" % ('qmu', self.qmu)
        s += "   %15s : %s\n" % ('autodiff', self.autodiff)
        s += "   %15s : %s\n" % ('smb', self.smb)
        #ENDDISP
        return s
    # }}}

    def VerboseToBinary(self):  # {{{
        #BEGINVERB2BIN
        binary = 0
        if self.mprocessor:
            binary = binary | 1
        if self.module:
            binary = binary | 2
        if self.solution:
            binary = binary | 4
        if self.solver:
            binary = binary | 8
        if self.convergence:
            binary = binary | 16
        if self.control:
            binary = binary | 32
        if self.qmu:
            binary = binary | 64
        if self.autodiff:
            binary = binary | 128
        if self.smb:
            binary = binary | 256
        #ENDVERB2BIN

        return binary
    # }}}

    def BinaryToVerbose(self, binary):  # {{{

        #BEGINBIN2VERB
        self.mprocessor = bool(binary & 1)
        self.module = bool(binary & 2)
        self.solution = bool(binary & 4)
        self.solver = bool(binary & 8)
        self.convergence = bool(binary & 16)
        self.control = bool(binary & 32)
        self.qmu = bool(binary & 64)
        self.autodiff = bool(binary & 128)
        self.smb = bool(binary & 256)
    #ENDBIN2VERB
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        return md
    # }}}

    def marshall(self, prefix, md, fid):  # {{{
        WriteData(fid, prefix, 'data', self.VerboseToBinary(), 'name', 'md.verbose', 'format', 'Integer')
    # }}}
