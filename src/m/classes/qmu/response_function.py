import numpy as np

from fielddisplay import fielddisplay
from MatlabFuncs import *
from pairoptions import pairoptions
from partition_npart import *
from rlev_write import *
#from rlist_write import *


class response_function(object):
    '''
    RESPONSE_FUNCTION class definition

        Usage:
            rf = response_function(
                'descriptor', descriptor,
                'response_levels', respl,
                'probability_levels', probl,
                'reliability_levels',rell,
                'general_reliability_levels', grell,
                'partition', partition
            )

            where rf is the response function object returned by the 
            constructor. All options except the descriptor are optional. A partition can be provided for scaled variables.

        Example:
            md.qmu.responses.maxvel = response_function(
                'descriptor', 'MaxVel',
                'response_levels',[0],
                'probl', [0.0001 0.001 0.01 0.25 0.5 0.75 0.99 0.999 0.9999]
                )
    '''
    def __init__(self):
        self.descriptor = ''
        self.respl      = []
        self.probl      = []
        self.rell       = []
        self.grell      = []
        self.partition  = []

    @staticmethod
    def response_function(*args):

        nargin = len(args)
        # create a default object
        if nargin == 0:
            return response_function()

        # copy the object or create the object from the input
        elif nargin == 1:
            if isinstance(args[0], response_function):
                rf = args[0]
            else:
                raise RuntimeError('Object ' + str(args[0]) + ' is a ' + str(type(args[0])) + ' class object, not "response_function".')
        else:
            # asizec = array_size(*args[0:min(nargin, 1)])
            # rf = [response_function() for i in range(asizec[0]) for j in range(asizec[1])]

            # for i in range(np.size(rf)):
            #     if (np.size(args[0]) > 1):
            #         rf[i].descriptor = args[0][i]
            #     else:
            #         rf[i].descriptor = str(args[0]) + string_dim(rf, i, 'vector')

            # if nargin >= 2:
            #     for i in range(np.size(rf)):
            #         rf[i].respl = args[1]

            # if nargin >= 3:
            #     for i in range(np.size(rf)):
            #         rf[i].probl = args[2]

            # if nargin >= 4:
            #     for i in range(np.size(rf)):
            #         rf[i].rell = args[3]

            # if nargin >= 5:
            #     for i in range(np.size(rf)):
            #         rf[i].grell = args[4]

            # if nargin > 5:
            #     print('WARNING: response_function:extra_arg: Extra arguments for object of class ' + str(type(rf)) + '.')

            rf = response_function()

            #recover options:
            options = pairoptions(*args)

            #initialize fields:
            rf.descriptor = options.getfieldvalue('descriptor')
            rf.respl = options.getfieldvalue('response_levels', [])
            rf.probl = options.getfieldvalue('probability_levels', [0.0001, 0.001, 0.01, 0.25, 0.5, 0.75, 0.99, 0.999, 0.9999])
            rf.rell = options.getfieldvalue('reliability_levels', [])
            rf.grell = options.getfieldvalue('general_reliability_levels', [])

            #if the response is scaled, a partition vector should have been supplied.
            if rf.isscaled():
                rf.partition = options.getfieldvalue('partition')
                npart = partition_npart(rf.partition)

        return [rf] # Always return a list, so we have something akin to a MATLAB single row matrix

    def __repr__(rf):  # {{{
        # display the object
        string = 'class "response_function" object = \n'
        string = "%s\n%s" % (string, fielddisplay(rf, 'descriptor', 'name tag'))
        string = "%s\n%s" % (string, fielddisplay(rf, 'respl', 'response levels'))
        string = "%s\n%s" % (string, fielddisplay(rf, 'probl', 'probability levels'))
        string = "%s\n%s" % (string, fielddisplay(rf, 'rell', 'reliability levels'))
        string = "%s\n%s" % (string, fielddisplay(rf, 'grell', 'general reliability levels'))

        if rf.partition != []:
            string = "%s\n%s" % (string, fielddisplay(rf, 'partition', 'partition vector defining where the response will be computed'))

        return string
    # }}}

    def __len__(self):  # {{{
        return max(len(self.respl), len(self.probl), len(self.rell), len(self.grell))
    # }}}

    @staticmethod
    def prop_desc(rf, dstr):  # {{{
        desc = ['' for i in range(np.size(rf))]
        for i in range(np.size(rf)):
            if rf[i].descriptor != '' or type(rf[i].descriptor) != str:
                desc[i] = str(rf[i].descriptor)
            elif dstr != '':
                desc[i] = str(dstr) + str(string_dim(rf, i, 'vector'))
            else:
                desc[i] = 'rf' + str(string_dim(rf, i, 'vector'))

        desc = allempty(desc)
        return desc
    # }}}

    @staticmethod
    def prop_stype(rf):  # {{{
        stype = []
        return stype
    # }}}

    @staticmethod
    def prop_scale(rf):  # {{{
        scale = []
        return scale
    # }}}

    @staticmethod
    def prop_weight(rf):  # {{{
        weight = []
        return weight
    # }}}

    @staticmethod
    def prop_lower(rf):  # {{{
        lower = []
        return lower
    # }}}

    @staticmethod
    def prop_upper(rf):  # {{{
        upper = []
        return upper
    # }}}

    @staticmethod
    def prop_target(rf):  # {{{
        target = []
        return target
    # }}}

    @staticmethod
    def prop_levels(rf):
        respl = empty_nd_list(np.size(rf))
        probl = empty_nd_list(np.size(rf))
        rell = empty_nd_list(np.size(rf))
        grell = empty_nd_list(np.size(rf))
        for i in range(np.size(rf)):
            respl[i] = rf[i].respl
            probl[i] = rf[i].probl
            rell[i] = rf[i].rell
            grell[i] = rf[i].grell
        respl = allempty(respl)
        probl = allempty(probl)
        rell = allempty(rell)
        grell = allempty(grell)
        return [respl, probl, rell, grell]

    #new methods:
    def isscaled(self):  # {{{
        if strncmpi(self.descriptor, 'scaled_', 7):
            return True
        else:
            return False
    # }}}

    @staticmethod
    def dakota_write(fidi, dresp, rdesc):
        #possible namespace pollution here
        from rlist_write import rlist_write
        # collect only the responses of the appropriate class
        #rf = [struc_class(vars(dresp)[i][j], 'response_function', 'rf') for i in fieldnames(dresp) for j in range(len(vars(dresp)[i]))]
        resp = deepcopy(dresp)
        fields = fieldnames(resp)
        for field in fields:
            if getattr(resp, field)[0].__class__.__name__ != 'response_function':
                delattr(resp, field)
        if len(resp) > 0:
            rdesc = rlist_write(fidi, 'response_function', 'rf', resp, rdesc)
        return rdesc

    @staticmethod
    def dakota_rlev_write(fidi, dresp, params):
        # collect only the responses of the appropriate class
        rf = [struc_class(vars(dresp)[i][j], 'response_function', 'rf') for i in fieldnames(dresp) for j in range(len(vars(dresp)[i]))]

        # write response levels
        rlev_write(fidi, rf, 'response_function', params)
