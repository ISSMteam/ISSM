
import numpy as np

from checkfield import checkfield
from fielddisplay import fielddisplay
import MatlabFuncs as m
from project3d import project3d
from WriteData import WriteData


class mmemasstransport(object):
    """
    MMEMASSTRANSPORT class definition

    Usage:
        masstransport=mmemasstransport();
    A little bit of info here. thickness is an ensemble set of element based time series (cell array 
    of matrice size nel x ntimes. 
    We supply indices into this ensemble: ex: [1 5 10] and a partition vector (values 0 to nindices-1). 
    This means, that when we resolve this MME, we'll have: 
    thickness = thickness{0} modified by: 
            - thickness{1} on elements such that their partition value is 0
            - thickness{5} on elements such that their partition value is 1
            - thickness{10} on elements such that their partition value is 2
    The goal is to allow to index into different Mmes. Here, partition value 0 might be for GlacierMIP2 
    elements, partition value 1 for ISMIP6 Antarctica and value 2 for ISMIP6 Greenland. 
    """

    def __init__(self,*args): # {{{
        self.thickness   = np.nan # multi model ensemble time series of ice thickness evolution
        self.ids         = [] # indices into the ensemble 
        self.partition   = [] # partition vector 
        self.requested_outputs = []
        # }}}

    def __repr__(self):  # {{{
        s = '   mmemasstransport:\n'
        s += '{}\n'.format(fielddisplay(self,'thickness','ensemble set of thickness changes'))
        s += '{}\n'.format(fielddisplay(self,'ids','indices  into the ensemble set'))
        s += '{}\n'.format(fielddisplay(self,'partition','partition vector'))
        s += '{}\n'.format(fielddisplay(self,'requested_outputs','additional outputs requested'))

        return s
        # }}}

    def extrude(self, md):  # {{{
        print('Nothing to be done with extrude in "mmemasstransport".')
        # }}}

    def setdefaultparameters(self): # {{{
        self.ids = []
        self.thickness = []
        self.partition = []
        #default output
        self.requested_outputs = ['default']
        # }}}
    
    def defaultoutputs(self,md): # {{{
        lists = ['DeltaIceThickness']
        return lists
        # }}}

    def checkconsistency(self, md, solution, analyses): # {{{

        #Early return
        if not 'MmemasstransportAnalysis' in analyses or \
            (solution == 'TransientSolution' and md.transient.ismmemasstransport==0):
            return
            
        nids=len(self.ids)
        nt=len(self.thickness)
        s = np.shape(self.ids)
        if s[0] < s[1]:
            raise Exception('mmemasstransport checkconsistency error message: ids should be a column vector');
        
        md = checkfield(md, 'field', self.partition, 'fieldname', 'partition', 'NaN', 1, 'Inf', 1, '>=', -1, '<=', nids-1, 'size', [md.mesh.numberofelements,1])
        md = checkfield(md, 'field', self.ids, 'fieldname', 'ids', 'NaN', 1, 'Inf', 1, '>=', 1, '<=', nt)
        for i in range(nt):
            md = checkfield(md, 'field', self.thickness[i], 'fieldname', 'self.thickness{' + str(i) + '}', 'NaN', 0, 'Inf', 1, 'timeseries', 1, 'timeserieslength', md.mesh.numberofelements+1)
             #don't check for NaN, as it's a legit value here for vertices that do not belong to the current Mme. 
        
        md = checkfield(md, 'fieldname', 'mmemasstransport.requested_outputs', 'stringrow', 1)
        # }}}

    def nummodels(self): # {{{
        return len(self.thickness)
        # }}}

    def marshall(self, prefix, md, fid): # {{{
        WriteData(fid,prefix,'object',self,'fieldname','ids','format','DoubleMat','mattype',3)
        WriteData(fid,prefix,'object',self,'fieldname','partition','format','DoubleMat','mattype',1)
        WriteData(fid,prefix,'object',self,'fieldname','thickness','format','MatArray','timeseries',1,'timeserieslength',md.mesh.numberofelements+1,'yts',md.constants.yts)

        #process requested outputs
        outputs = self.requested_outputs
        indices = [i for i, x in enumerate(outputs) if x == 'default']
        if len(indices) > 0:
            outputscopy = outputs[0:max(0, indices[0] - 1)] + self.defaultoutputs(md) + outputs[indices[0] + 1:]
            outputs = outputscopy
    
        WriteData(fid,prefix,'data',outputs,'name','md.mmemasstransport.requested_outputs','format','StringArray')
        # }}}


