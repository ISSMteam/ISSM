#module imports {{{
import numpy as np
import copy
import sys
from pairoptions import pairoptions
from mesh2d import mesh2d
from mesh3dprisms import mesh3dprisms
from mask import mask
from geometry import geometry
from constants import constants
from SMBforcing import SMBforcing
from SMBpdd import SMBpdd
from SMBd18opdd import SMBd18opdd
from SMBgradients import SMBgradients
from SMBcomponents import SMBcomponents
from SMBmeltcomponents import SMBmeltcomponents
from basalforcings import basalforcings
from linearbasalforcings import linearbasalforcings
from matice import matice
from levelset import levelset
from calving import calving
from love import love
from calvinglevermann import calvinglevermann
#from calvingpi import calvingpi
from frontalforcings import frontalforcings
from damage import damage
from friction import friction
from flowequation import flowequation
from timestepping import timestepping
from timesteppingadaptive import timesteppingadaptive
from initialization import initialization
from rifts import rifts
from solidearth import solidearth
from dsl import dsl
from debug import debug
from verbose import verbose
from issmsettings import issmsettings
from toolkits import toolkits
from generic import generic
from balancethickness import balancethickness
from stressbalance import stressbalance
from groundingline import groundingline
from hydrologyshreve import hydrologyshreve
from hydrologydc import hydrologydc
from hydrologyglads import hydrologyglads
from hydrologypism import hydrologypism
from hydrologyshakti import hydrologyshakti
from debris import debris
from masstransport import masstransport
from thermal import thermal
from steadystate import steadystate
from transient import transient
from esa import esa
from autodiff import autodiff
from inversion import inversion
from outputdefinition import outputdefinition
from qmu import qmu
from amr import amr
from results import results
from radaroverlay import radaroverlay
from miscellaneous import miscellaneous
from private import private
from mumpsoptions import mumpsoptions
from iluasmoptions import iluasmoptions
from project3d import project3d
from project2d import project2d
from FlagElements import FlagElements
from NodeConnectivity import NodeConnectivity
from ElementConnectivity import ElementConnectivity
from contourenvelope import contourenvelope
from DepthAverage import DepthAverage
from sampling import sampling
from stochasticforcing import stochasticforcing
# }}}


class model(object):
    """model class definition

    Usage:
        md = model()
    """

    def __init__(self, *args):  #{{{
        self.mesh = None
        self.mask = None

        self.geometry = None
        self.constants = None
        self.smb = None
        self.basalforcings = None
        self.materials = None
        self.damage = None
        self.friction = None
        self.flowequation = None
        self.timestepping = None
        self.initialization = None
        self.rifts = None
        self.dsl = None
        self.solidearth = None
        self.debug = None
        self.verbose = None
        self.settings = None
        self.toolkits = None
        self.cluster = None
        self.balancethickness = None
        self.stressbalance = None
        self.groundingline = None
        self.hydrology = None
        self.debris = None
        self.masstransport = None
        self.thermal = None
        self.steadystate = None
        self.transient = None
        self.levelset = None
        self.calving = None
        self.frontalforcings = None
        self.love = None
        self.esa = None
        self.sampling = None
        self.autodiff = None
        self.inversion = None
        self.qmu = None
        self.amr = None
        self.results = None
        self.outputdefinition = None
        self.radaroverlay = None
        self.miscellaneous = None
        self.private = None
        self.stochasticforcing = None

        if len(args) == 0:
            self.setdefaultparameters('earth')
        else:
            options = pairoptions(*args)
            planet = options.getfieldvalue('planet', 'earth')
            self.setdefaultparameters(planet)
    # }}}

    def __repr__(obj):  #{{{
        # TODO:
        # - Convert all formatting to calls to <string>.format (see any
        #   already converted <class>.__repr__ method for examples)
        #
        s = '%19s: %-23s -- %s' % ('mesh', '[%s %s]' % ('1x1', obj.mesh.__class__.__name__), 'mesh properties')
        s = '%s\n%s' % (s, '%19s: %-23s -- %s' % ('mask', '[%s %s]' % ('1x1', obj.mask.__class__.__name__), 'defines grounded and floating elements'))
        s = '%s\n%s' % (s, '%19s: %-23s -- %s' % ('geometry', '[%s %s]' % ('1x1', obj.geometry.__class__.__name__), 'surface elevation, bedrock topography, ice thickness, ...'))
        s = '%s\n%s' % (s, '%19s: %-23s -- %s' % ('constants', '[%s %s]' % ('1x1', obj.constants.__class__.__name__), 'physical constants'))
        s = '%s\n%s' % (s, '%19s: %-23s -- %s' % ('smb', '[%s %s]' % ('1x1', obj.smb.__class__.__name__), 'surface mass balance'))
        s = '%s\n%s' % (s, '%19s: %-23s -- %s' % ('basalforcings', '[%s %s]' % ('1x1', obj.basalforcings.__class__.__name__), 'bed forcings'))
        s = '%s\n%s' % (s, '%19s: %-23s -- %s' % ('materials', '[%s %s]' % ('1x1', obj.materials.__class__.__name__), 'material properties'))
        s = '%s\n%s' % (s, '%19s: %-23s -- %s' % ('damage', '[%s %s]' % ('1x1', obj.damage.__class__.__name__), 'damage propagation laws'))
        s = '%s\n%s' % (s, '%19s: %-23s -- %s' % ('friction', '[%s %s]' % ('1x1', obj.friction.__class__.__name__), 'basal friction / drag properties'))
        s = '%s\n%s' % (s, '%19s: %-23s -- %s' % ('flowequation', '[%s %s]' % ('1x1', obj.flowequation.__class__.__name__), 'flow equations'))
        s = '%s\n%s' % (s, '%19s: %-23s -- %s' % ('timestepping', '[%s %s]' % ('1x1', obj.timestepping.__class__.__name__), 'time stepping for transient models'))
        s = '%s\n%s' % (s, '%19s: %-23s -- %s' % ('initialization', '[%s %s]' % ('1x1', obj.initialization.__class__.__name__), 'initial guess / state'))
        s = '%s\n%s' % (s, '%19s: %-23s -- %s' % ('rifts', '[%s %s]' % ('1x1', obj.rifts.__class__.__name__), 'rifts properties'))
        s = '%s\n%s' % (s, '%19s: %-23s -- %s' % ('solidearth', '[%s %s]' % ('1x1', obj.solidearth.__class__.__name__), 'solidearth inputs and settings'))
        s = '%s\n%s' % (s, '%19s: %-23s -- %s' % ('dsl', '[%s %s]' % ('1x1', obj.dsl.__class__.__name__), 'dynamic sea level'))
        s = '%s\n%s' % (s, '%19s: %-23s -- %s' % ('debug', '[%s %s]' % ('1x1', obj.debug.__class__.__name__), 'debugging tools (valgrind, gprof)'))
        s = '%s\n%s' % (s, '%19s: %-23s -- %s' % ('verbose', '[%s %s]' % ('1x1', obj.verbose.__class__.__name__), 'verbosity level in solve'))
        s = '%s\n%s' % (s, '%19s: %-23s -- %s' % ('settings', '[%s %s]' % ('1x1', obj.settings.__class__.__name__), 'settings properties'))
        s = '%s\n%s' % (s, '%19s: %-23s -- %s' % ('toolkits', '[%s %s]' % ('1x1', obj.toolkits.__class__.__name__), 'PETSc options for each solution'))
        s = '%s\n%s' % (s, '%19s: %-23s -- %s' % ('cluster', '[%s %s]' % ('1x1', obj.cluster.__class__.__name__), 'cluster parameters (number of CPUs...)'))
        s = '%s\n%s' % (s, '%19s: %-23s -- %s' % ('balancethickness', '[%s %s]' % ('1x1', obj.balancethickness.__class__.__name__), 'parameters for balancethickness solution'))
        s = '%s\n%s' % (s, '%19s: %-23s -- %s' % ('stressbalance', '[%s %s]' % ('1x1', obj.stressbalance.__class__.__name__), 'parameters for stressbalance solution'))
        s = '%s\n%s' % (s, '%19s: %-23s -- %s' % ('groundingline', '[%s %s]' % ('1x1', obj.groundingline.__class__.__name__), 'parameters for groundingline solution'))
        s = '%s\n%s' % (s, '%19s: %-23s -- %s' % ('hydrology', '[%s %s]' % ('1x1', obj.hydrology.__class__.__name__), 'parameters for hydrology solution'))
        s = '%s\n%s' % (s, '%19s: %-23s -- %s' % ('masstransport', '[%s %s]' % ('1x1', obj.masstransport.__class__.__name__), 'parameters for masstransport solution'))
        s = '%s\n%s' % (s, '%19s: %-23s -- %s' % ('thermal', '[%s %s]' % ('1x1', obj.thermal.__class__.__name__), 'parameters for thermal solution'))
        s = '%s\n%s' % (s, '%19s: %-23s -- %s' % ('steadystate', '[%s %s]' % ('1x1', obj.steadystate.__class__.__name__), 'parameters for steadystate solution'))
        s = '%s\n%s' % (s, '%19s: %-23s -- %s' % ('transient', '[%s %s]' % ('1x1', obj.transient.__class__.__name__), 'parameters for transient solution'))
        s = '%s\n%s' % (s, '%19s: %-23s -- %s' % ('levelset', '[%s %s]' % ('1x1', obj.levelset.__class__.__name__), 'parameters for moving boundaries (level-set method)'))
        s = '%s\n%s' % (s, '%19s: %-23s -- %s' % ('calving', '[%s %s]' % ('1x1', obj.calving.__class__.__name__), 'parameters for calving'))
        s = '%s\n%s' % (s, '%19s: %-23s -- %s' % ('frontalforcings', '[%s %s]' % ('1x1', obj.frontalforcings.__class__.__name__), 'parameters for frontalforcings'))
        s = '%s\n%s' % (s, '%19s: %-23s -- %s' % ('esa', '[%s %s]' % ('1x1', obj.esa.__class__.__name__), 'parameters for elastic adjustment solution'))
        s = '%s\n%s' % (s, '%19s: %-23s -- %s' % ('sampling', '[%s %s]' % ('1x1', obj.sampling.__class__.__name__), 'parameters for stochastic sampler'))
        s = '%s\n%s' % (s, '%19s: %-23s -- %s' % ('love', '[%s %s]' % ('1x1', obj.love.__class__.__name__), 'parameters for love solution'))
        s = '%s\n%s' % (s, '%19s: %-23s -- %s' % ('autodiff', '[%s %s]' % ('1x1', obj.autodiff.__class__.__name__), 'automatic differentiation parameters'))
        s = '%s\n%s' % (s, '%19s: %-23s -- %s' % ('inversion', '[%s %s]' % ('1x1', obj.inversion.__class__.__name__), 'parameters for inverse methods'))
        s = '%s\n%s' % (s, '%19s: %-23s -- %s' % ('qmu', '[%s %s]' % ('1x1', obj.qmu.__class__.__name__), 'Dakota properties'))
        s = '%s\n%s' % (s, '%19s: %-23s -- %s' % ('amr', '[%s %s]' % ('1x1', obj.amr.__class__.__name__), 'adaptive mesh refinement properties'))
        s = '%s\n%s' % (s, '%19s: %-23s -- %s' % ('outputdefinition', '[%s %s]' % ('1x1', obj.outputdefinition.__class__.__name__), 'output definition'))
        s = '%s\n%s' % (s, '%19s: %-23s -- %s' % ('results', '[%s %s]' % ('1x1', obj.results.__class__.__name__), 'model results'))
        s = '%s\n%s' % (s, '%19s: %-23s -- %s' % ('radaroverlay', '[%s %s]' % ('1x1', obj.radaroverlay.__class__.__name__), 'radar image for plot overlay'))
        s = '%s\n%s' % (s, '%19s: %-23s -- %s' % ('miscellaneous', '[%s %s]' % ('1x1', obj.miscellaneous.__class__.__name__), 'miscellaneous fields'))
        s = '%s\n%s' % (s, '%19s: %-23s -- %s' % ('stochasticforcing', '[%s %s]' % ('1x1', obj.stochasticforcing.__class__.__name__), 'stochasticity applied to model forcings'))
        return s
    # }}}

    def properties(self):  #{{{
        # ordered list of properties since vars(self) is random
        return [
            'mesh',
            'mask',
            'geometry',
            'constants',
            'smb',
            'basalforcings',
            'materials',
            'damage',
            'friction',
            'flowequation',
            'timestepping',
            'initialization',
            'rifts',
            'dsl',
            'solidearth',
            'debug',
            'verbose',
            'settings',
            'toolkits',
            'cluster',
            'balancethickness',
            'stressbalance',
            'groundingline',
            'hydrology',
            'debris',
            'masstransport',
            'thermal',
            'steadystate',
            'transient',
            'levelset',
            'calving',
            'frontalforcings',
            'love',
            'esa',
            'sampling',
            'autodiff',
            'inversion',
            'qmu',
            'amr',
            'results',
            'outputdefinition',
            'radaroverlay',
            'miscellaneous',
            'private',
            'stochasticforcing'
        ]
    # }}}

    def setdefaultparameters(self, planet):  #{{{
        self.mesh = mesh2d()
        self.mask = mask()
        self.constants = constants()
        self.geometry = geometry()
        self.initialization = initialization()
        self.smb = SMBforcing()
        self.basalforcings = basalforcings()
        self.friction = friction()
        self.rifts = rifts()
        self.solidearth = solidearth(planet)
        self.dsl = dsl()
        self.timestepping = timestepping()
        self.groundingline = groundingline()
        self.materials = matice()
        self.damage = damage()
        self.flowequation = flowequation()
        self.debug = debug()
        self.verbose = verbose()
        self.settings = issmsettings()
        self.toolkits = toolkits()
        self.cluster = generic()
        self.balancethickness = balancethickness()
        self.stressbalance = stressbalance()
        self.hydrology = hydrologyshreve()
        self.debris = debris()
        self.masstransport = masstransport()
        self.thermal = thermal()
        self.steadystate = steadystate()
        self.transient = transient()
        self.levelset = levelset()
        self.calving = calving()
        self.frontalforcings = frontalforcings()
        self.love = love()
        self.esa = esa()
        self.sampling = sampling()
        self.autodiff = autodiff()
        self.inversion = inversion()
        self.qmu = qmu()
        self.amr = amr()
        self.radaroverlay = radaroverlay()
        self.results = results()
        self.outputdefinition = outputdefinition()
        self.miscellaneous = miscellaneous()
        self.private = private()
        self.stochasticforcing = stochasticforcing()
    # }}}

    def checkmessage(self, string):  #{{{
        print("model not consistent: {}".format(string))
        self.private.isconsistent = False
        return self
    # }}}
    #@staticmethod

    def extract(self, area):  #{{{
        """EXTRACT - extract a model according to an Argus contour or flag list

        This routine extracts a submodel from a bigger model with respect to a given contour
        md must be followed by the corresponding exp file or flags list
        It can either be a domain file (argus type, .exp extension), or an array of element flags.
        If user wants every element outside the domain to be
        extract2d, add '~' to the name of the domain file (ex: '~HO.exp')
        an empty string '' will be considered as an empty domain
        a string 'all' will be considered as the entire domain

        Usage:
            md2 = extract(md, area)

        Examples:
            md2 = extract(md, 'Domain.exp')

        See also: EXTRUDE, COLLAPSE
        """

        #copy model
        md1 = copy.deepcopy(self)

        #get elements that are inside area
        flag_elem = FlagElements(md1, area)
        if not np.any(flag_elem):
            raise RuntimeError("extracted model is empty")

        #kick out all elements with 3 dirichlets
        spc_elem = np.nonzero(np.logical_not(flag_elem))[0]
        spc_node = np.unique(md1.mesh.elements[spc_elem, :]) - 1
        flag = np.ones(md1.mesh.numberofvertices)
        flag[spc_node] = 0
        pos = np.nonzero(np.logical_not(np.sum(flag[md1.mesh.elements - 1], axis=1)))[0]
        flag_elem[pos] = 0

        #extracted elements and nodes lists
        pos_elem = np.nonzero(flag_elem)[0]
        pos_node = np.unique(md1.mesh.elements[pos_elem, :]) - 1

        #keep track of some fields
        numberofvertices1 = md1.mesh.numberofvertices
        numberofelements1 = md1.mesh.numberofelements
        numberofvertices2 = np.size(pos_node)
        numberofelements2 = np.size(pos_elem)
        flag_node = np.zeros(numberofvertices1)
        flag_node[pos_node] = 1

        #Create Pelem and Pnode (transform old nodes in new nodes and same thing for the elements)
        Pelem = np.zeros(numberofelements1, int)
        Pelem[pos_elem] = np.arange(1, numberofelements2 + 1)
        Pnode = np.zeros(numberofvertices1, int)
        Pnode[pos_node] = np.arange(1, numberofvertices2 + 1)

        #renumber the elements (some node won't exist anymore)
        elements_1 = copy.deepcopy(md1.mesh.elements)
        elements_2 = elements_1[pos_elem, :]
        elements_2[:, 0] = Pnode[elements_2[:, 0] - 1]
        elements_2[:, 1] = Pnode[elements_2[:, 1] - 1]
        elements_2[:, 2] = Pnode[elements_2[:, 2] - 1]
        if md1.mesh.__class__.__name__ == 'mesh3dprisms':
            elements_2[:, 3] = Pnode[elements_2[:, 3] - 1]
            elements_2[:, 4] = Pnode[elements_2[:, 4] - 1]
            elements_2[:, 5] = Pnode[elements_2[:, 5] - 1]

        #OK, now create the new model!

        #take every field from model
        md2 = copy.deepcopy(md1)

        #automatically modify fields

        #loop over model fields
        model_fields = vars(md1)
        for fieldi in model_fields:
            #get field
            field = getattr(md1, fieldi)
            fieldsize = np.shape(field)
            if hasattr(field, '__dict__') and fieldi not in ['results']:  #recursive call
                object_fields = vars(field)
                for fieldj in object_fields:
                    #get field
                    field = getattr(getattr(md1, fieldi), fieldj)
                    fieldsize = np.shape(field)
                    if len(fieldsize):
                        #size = number of nodes * n
                        if fieldsize[0] == numberofvertices1:
                            setattr(getattr(md2, fieldi), fieldj, field[pos_node])
                        elif fieldsize[0] == numberofvertices1 + 1:
                            setattr(getattr(md2, fieldi), fieldj, np.vstack((field[pos_node], field[-1, :])))
                        #size = number of elements * n
                        elif fieldsize[0] == numberofelements1:
                            setattr(getattr(md2, fieldi), fieldj, field[pos_elem])
            else:
                if len(fieldsize):
                    #size = number of nodes * n
                    if fieldsize[0] == numberofvertices1:
                        setattr(md2, fieldi, field[pos_node])
                    elif fieldsize[0] == numberofvertices1 + 1:
                        setattr(md2, fieldi, np.hstack((field[pos_node], field[-1, :])))
                    #size = number of elements * n
                    elif fieldsize[0] == numberofelements1:
                        setattr(md2, fieldi, field[pos_elem])

        #modify some specific fields
        #Mesh
        md2.mesh.numberofelements = numberofelements2
        md2.mesh.numberofvertices = numberofvertices2
        md2.mesh.elements = elements_2

        #mesh.uppervertex mesh.lowervertex
        if md1.mesh.__class__.__name__ == 'mesh3dprisms':
            md2.mesh.uppervertex = md1.mesh.uppervertex[pos_node]
            pos = np.where(~np.isnan(md2.mesh.uppervertex))[0]
            md2.mesh.uppervertex[pos] = Pnode[md2.mesh.uppervertex[pos].astype(int) - 1]

            md2.mesh.lowervertex = md1.mesh.lowervertex[pos_node]
            pos = np.where(~np.isnan(md2.mesh.lowervertex))[0]
            md2.mesh.lowervertex[pos] = Pnode[md2.mesh.lowervertex[pos].astype(int) - 1]

            md2.mesh.upperelements = md1.mesh.upperelements[pos_elem]
            pos = np.where(~np.isnan(md2.mesh.upperelements))[0]
            md2.mesh.upperelements[pos] = Pelem[md2.mesh.upperelements[pos].astype(int) - 1]

            md2.mesh.lowerelements = md1.mesh.lowerelements[pos_elem]
            pos = np.where(~np.isnan(md2.mesh.lowerelements))[0]
            md2.mesh.lowerelements[pos] = Pelem[md2.mesh.lowerelements[pos].astype(int) - 1]

        #Initial 2d mesh
        if md1.mesh.__class__.__name__ == 'mesh3dprisms':
            flag_elem_2d = flag_elem[np.arange(0, md1.mesh.numberofelements2d)]
            pos_elem_2d = np.nonzero(flag_elem_2d)[0]
            flag_node_2d = flag_node[np.arange(0, md1.mesh.numberofvertices2d)]
            pos_node_2d = np.nonzero(flag_node_2d)[0]

            md2.mesh.numberofelements2d = np.size(pos_elem_2d)
            md2.mesh.numberofvertices2d = np.size(pos_node_2d)
            md2.mesh.elements2d = md1.mesh.elements2d[pos_elem_2d, :]
            md2.mesh.elements2d[:, 0] = Pnode[md2.mesh.elements2d[:, 0] - 1]
            md2.mesh.elements2d[:, 1] = Pnode[md2.mesh.elements2d[:, 1] - 1]
            md2.mesh.elements2d[:, 2] = Pnode[md2.mesh.elements2d[:, 2] - 1]

            md2.mesh.x2d = md1.mesh.x[pos_node_2d]
            md2.mesh.y2d = md1.mesh.y[pos_node_2d]

        #Edges
        if md1.mesh.domaintype() == '2Dhorizontal':
            if np.ndim(md2.mesh.edges) > 1 and np.size(md2.mesh.edges, axis=1) > 1:  #do not use ~isnan because there are some np.nans...
                #renumber first two columns
                pos = np.nonzero(md2.mesh.edges[:, 3] != -1)[0]
                md2.mesh.edges[:, 0] = Pnode[md2.mesh.edges[:, 0] - 1]
                md2.mesh.edges[:, 1] = Pnode[md2.mesh.edges[:, 1] - 1]
                md2.mesh.edges[:, 2] = Pelem[md2.mesh.edges[:, 2] - 1]
                md2.mesh.edges[pos, 3] = Pelem[md2.mesh.edges[pos, 3] - 1]
                #remove edges when the 2 vertices are not in the domain.
                md2.mesh.edges = md2.mesh.edges[np.nonzero(np.logical_and(md2.mesh.edges[:, 0], md2.mesh.edges[:, 1]))[0], :]
                #Replace all zeros by - 1 in the last two columns
                pos = np.nonzero(md2.mesh.edges[:, 2] == 0)[0]
                md2.mesh.edges[pos, 2] = -1
                pos = np.nonzero(md2.mesh.edges[:, 3] == 0)[0]
                md2.mesh.edges[pos, 3] = -1
                #Invert - 1 on the third column with last column (Also invert first two columns!!)
                pos = np.nonzero(md2.mesh.edges[:, 2] == -1)[0]
                md2.mesh.edges[pos, 2] = md2.mesh.edges[pos, 3]
                md2.mesh.edges[pos, 3] = -1
                values = md2.mesh.edges[pos, 1]
                md2.mesh.edges[pos, 1] = md2.mesh.edges[pos, 0]
                md2.mesh.edges[pos, 0] = values
                #Finally remove edges that do not belong to any element
                pos = np.nonzero(np.logical_and(md2.mesh.edges[:, 1] == -1, md2.mesh.edges[:, 2] == -1))[0]
                md2.mesh.edges = np.delete(md2.mesh.edges, pos, axis=0)

        #Penalties
        if np.any(np.logical_not(np.isnan(md2.stressbalance.vertex_pairing))):
            for i in range(np.size(md1.stressbalance.vertex_pairing, axis=0)):
                md2.stressbalance.vertex_pairing[i, :] = Pnode[md1.stressbalance.vertex_pairing[i, :]]
            md2.stressbalance.vertex_pairing = md2.stressbalance.vertex_pairing[np.nonzero(md2.stressbalance.vertex_pairing[:, 0])[0], :]
        if np.any(np.logical_not(np.isnan(md2.masstransport.vertex_pairing))):
            for i in range(np.size(md1.masstransport.vertex_pairing, axis=0)):
                md2.masstransport.vertex_pairing[i, :] = Pnode[md1.masstransport.vertex_pairing[i, :]]
            md2.masstransport.vertex_pairing = md2.masstransport.vertex_pairing[np.nonzero(md2.masstransport.vertex_pairing[:, 0])[0], :]

        #recreate segments
        if md1.mesh.__class__.__name__ == 'mesh2d':
            md2.mesh.vertexconnectivity = NodeConnectivity(md2.mesh.elements, md2.mesh.numberofvertices)
            md2.mesh.elementconnectivity = ElementConnectivity(md2.mesh.elements, md2.mesh.vertexconnectivity)
            md2.mesh.segments = contourenvelope(md2.mesh)
            md2.mesh.vertexonboundary = np.zeros(numberofvertices2, bool)
            md2.mesh.vertexonboundary[md2.mesh.segments[:, 0:2] - 1] = True
        else:
            #First do the connectivity for the contourenvelope in 2d
            md2.mesh.vertexconnectivity = NodeConnectivity(md2.mesh.elements2d, md2.mesh.numberofvertices2d)
            md2.mesh.elementconnectivity = ElementConnectivity(md2.mesh.elements2d, md2.mesh.vertexconnectivity)
            segments = contourenvelope(md2.mesh)
            md2.mesh.vertexonboundary = np.zeros(int(numberofvertices2 / md2.mesh.numberoflayers), bool)
            md2.mesh.vertexonboundary[segments[:, 0:2] - 1] = True
            md2.mesh.vertexonboundary = np.tile(md2.mesh.vertexonboundary, md2.mesh.numberoflayers)
            #Then do it for 3d as usual
            md2.mesh.vertexconnectivity = NodeConnectivity(md2.mesh.elements, md2.mesh.numberofvertices)
            md2.mesh.elementconnectivity = ElementConnectivity(md2.mesh.elements, md2.mesh.vertexconnectivity)

        #Boundary conditions: Dirichlets on new boundary
        #Catch the elements that have not been extracted
        orphans_elem = np.nonzero(np.logical_not(flag_elem))[0]
        orphans_node = np.unique(md1.mesh.elements[orphans_elem, :]) - 1
        #Figure out which node are on the boundary between md2 and md1
        nodestoflag1 = np.intersect1d(orphans_node, pos_node)
        nodestoflag2 = Pnode[nodestoflag1].astype(int) - 1
        if np.size(md1.stressbalance.spcvx) > 1 and np.size(md1.stressbalance.spcvy) > 2 and np.size(md1.stressbalance.spcvz) > 2:
            if np.size(md1.inversion.vx_obs) > 1 and np.size(md1.inversion.vy_obs) > 1:
                md2.stressbalance.spcvx[nodestoflag2] = md2.inversion.vx_obs[nodestoflag2]
                md2.stressbalance.spcvy[nodestoflag2] = md2.inversion.vy_obs[nodestoflag2]
            else:
                md2.stressbalance.spcvx[nodestoflag2] = np.nan
                md2.stressbalance.spcvy[nodestoflag2] = np.nan
                print("\n!! extract warning: spc values should be checked !!\n\n")
            #put 0 for vz
            md2.stressbalance.spcvz[nodestoflag2] = 0
        if np.any(np.logical_not(np.isnan(md1.thermal.spctemperature))):
            md2.thermal.spctemperature[nodestoflag2] = 1

        #Results fields
        if md1.results:
            md2.results = results()
            for solutionfield, field in list(md1.results.__dict__.items()):
                if isinstance(field, list):
                    setattr(md2.results, solutionfield, [])
                    #get time step
                    for i, fieldi in enumerate(field):
                        if isinstance(fieldi, results) and fieldi:
                            getattr(md2.results, solutionfield).append(results())
                            fieldr = getattr(md2.results, solutionfield)[i]
                            #get subfields
                            for solutionsubfield, subfield in list(fieldi.__dict__.items()):
                                if np.size(subfield) == numberofvertices1:
                                    setattr(fieldr, solutionsubfield, subfield[pos_node])
                                elif np.size(subfield) == numberofelements1:
                                    setattr(fieldr, solutionsubfield, subfield[pos_elem])
                                else:
                                    setattr(fieldr, solutionsubfield, subfield)
                        else:
                            getattr(md2.results, solutionfield).append(None)
                elif isinstance(field, results):
                    setattr(md2.results, solutionfield, results())
                    if isinstance(field, results) and field:
                        fieldr = getattr(md2.results, solutionfield)
                        #get subfields
                        for solutionsubfield, subfield in list(field.__dict__.items()):
                            if np.size(subfield) == numberofvertices1:
                                setattr(fieldr, solutionsubfield, subfield[pos_node])
                            elif np.size(subfield) == numberofelements1:
                                setattr(fieldr, solutionsubfield, subfield[pos_elem])
                            else:
                                setattr(fieldr, solutionsubfield, subfield)

        #OutputDefinitions fields
        if md1.outputdefinition.definitions:
            for solutionfield, field in list(md1.outputdefinition.__dict__.items()):
                if isinstance(field, list):
                    #get each definition
                    for i, fieldi in enumerate(field):
                        if fieldi:
                            fieldr = getattr(md2.outputdefinition, solutionfield)[i]
                            #get subfields
                            for solutionsubfield, subfield in list(fieldi.__dict__.items()):
                                if np.size(subfield) == numberofvertices1:
                                    setattr(fieldr, solutionsubfield, subfield[pos_node])
                                elif np.size(subfield) == numberofelements1:
                                    setattr(fieldr, solutionsubfield, subfield[pos_elem])
                                else:
                                    setattr(fieldr, solutionsubfield, subfield)

        #Keep track of pos_node and pos_elem
        md2.mesh.extractedvertices = pos_node + 1
        md2.mesh.extractedelements = pos_elem + 1

        return md2
    # }}}

    def extrude(md, *args):  #{{{
        """EXTRUDE - vertically extrude a 2d mesh

        vertically extrude a 2d mesh and create corresponding 3d mesh.
        The vertical distribution can:
        - follow a polynomial law
        - follow two polynomial laws, one for the lower part and one for the upper part of the mesh
        - be discribed by a list of coefficients (between 0 and 1)


        Usage:
            md = extrude(md, numlayers, extrusionexponent)
            md = extrude(md, numlayers, lowerexponent, upperexponent)
            md = extrude(md, listofcoefficients)

        Example:
            md = extrude(md, 15, 1.3)
            md = extrude(md, 15, 1.3, 1.2)
            md = extrude(md, [0 0.2 0.5 0.7 0.9 0.95 1])

        See also: MODELEXTRACT, COLLAPSE
        """

        #some checks on list of arguments
        if len(args) > 3 or len(args) < 1:
            raise RuntimeError('extrude error message')

        #Extrude the mesh
        if len(args) == 1:  #list of coefficients
            clist = args[0]
            if any(clist < 0) or any(clist > 1):
                raise TypeError('extrusioncoefficients must be between 0 and 1')
            clist.extend([0., 1.])
            clist.sort()
            extrusionlist = list(set(clist))
            numlayers = len(extrusionlist)

        elif len(args) == 2:  #one polynomial law
            if args[1] <= 0:
                raise TypeError('extrusionexponent must be >= 0')
            numlayers = args[0]
            extrusionlist = (np.arange(0., float(numlayers - 1) + 1., 1.) / float(numlayers - 1))**args[1]

        elif len(args) == 3:  #two polynomial laws
            numlayers = args[0]
            lowerexp = args[1]
            upperexp = args[2]

            if args[1] <= 0 or args[2] <= 0:
                raise TypeError('lower and upper extrusionexponents must be >= 0')

            lowerextrusionlist = (np.arange(0., 1. + 2. / float(numlayers - 1), 2. / float(numlayers - 1)))**lowerexp / 2.
            upperextrusionlist = (np.arange(0., 1. + 2. / float(numlayers - 1), 2. / float(numlayers - 1)))**upperexp / 2.
            extrusionlist = np.unique(np.concatenate((lowerextrusionlist, 1. - upperextrusionlist)))

        if numlayers < 2:
            raise TypeError('number of layers should be at least 2')
        if md.mesh.__class__.__name__ == 'mesh3dprisms':
            raise TypeError('Cannot extrude a 3d mesh (extrude cannot be called more than once)')

        #Initialize with 2d mesh
        mesh2d = md.mesh
        md.mesh = mesh3dprisms()
        md.mesh.x = mesh2d.x
        md.mesh.y = mesh2d.y
        md.mesh.elements = mesh2d.elements
        md.mesh.numberofelements = mesh2d.numberofelements
        md.mesh.numberofvertices = mesh2d.numberofvertices

        md.mesh.lat = mesh2d.lat
        md.mesh.long = mesh2d.long
        md.mesh.epsg = mesh2d.epsg
        md.mesh.scale_factor = mesh2d.scale_factor

        md.mesh.vertexonboundary = mesh2d.vertexonboundary
        md.mesh.vertexconnectivity = mesh2d.vertexconnectivity
        md.mesh.elementconnectivity = mesh2d.elementconnectivity
        md.mesh.average_vertex_connectivity = mesh2d.average_vertex_connectivity

        md.mesh.extractedvertices = mesh2d.extractedvertices
        md.mesh.extractedelements = mesh2d.extractedelements

        x3d = np.empty((0))
        y3d = np.empty((0))
        z3d = np.empty((0))  #the lower node is on the bed
        thickness3d = md.geometry.thickness  #thickness and bed for these nodes
        bed3d = md.geometry.base

        #Create the new layers
        for i in range(numlayers):
            x3d = np.concatenate((x3d, md.mesh.x))
            y3d = np.concatenate((y3d, md.mesh.y))
            #nodes are distributed between bed and surface accordingly to the given exponent
            z3d = np.concatenate((z3d, (bed3d + thickness3d * extrusionlist[i]).reshape(-1)))
        number_nodes3d = np.size(x3d)  #number of 3d nodes for the non extruded part of the mesh

        #Extrude elements
        elements3d = np.empty((0, 6), int)
        for i in range(numlayers - 1):
            elements3d = np.vstack((elements3d, np.hstack((md.mesh.elements + i * md.mesh.numberofvertices,
                                                           md.mesh.elements + (i + 1) * md.mesh.numberofvertices))))  #Create the elements of the 3d mesh for the non extruded part
        number_el3d = np.size(elements3d, axis=0)  #number of 3d nodes for the non extruded part of the mesh

        #Keep a trace of lower and upper nodes
        lowervertex = np.nan * np.ones(number_nodes3d, int)
        uppervertex = np.nan * np.ones(number_nodes3d, int)
        lowervertex[md.mesh.numberofvertices:] = np.arange(1, (numlayers - 1) * md.mesh.numberofvertices + 1)
        uppervertex[:(numlayers - 1) * md.mesh.numberofvertices] = np.arange(md.mesh.numberofvertices + 1, number_nodes3d + 1)
        md.mesh.lowervertex = lowervertex
        md.mesh.uppervertex = uppervertex

        #same for lower and upper elements
        lowerelements = np.nan * np.ones(number_el3d, int)
        upperelements = np.nan * np.ones(number_el3d, int)
        lowerelements[md.mesh.numberofelements:] = np.arange(1, (numlayers - 2) * md.mesh.numberofelements + 1)
        upperelements[:(numlayers - 2) * md.mesh.numberofelements] = np.arange(md.mesh.numberofelements + 1, (numlayers - 1) * md.mesh.numberofelements + 1)
        md.mesh.lowerelements = lowerelements
        md.mesh.upperelements = upperelements

        #Save old mesh
        md.mesh.x2d = md.mesh.x
        md.mesh.y2d = md.mesh.y
        md.mesh.elements2d = md.mesh.elements
        md.mesh.numberofelements2d = md.mesh.numberofelements
        md.mesh.numberofvertices2d = md.mesh.numberofvertices

        #Build global 3d mesh
        md.mesh.elements = elements3d
        md.mesh.x = x3d
        md.mesh.y = y3d
        md.mesh.z = z3d
        md.mesh.numberofelements = number_el3d
        md.mesh.numberofvertices = number_nodes3d
        md.mesh.numberoflayers = numlayers

        #Ok, now deal with the other fields from the 2d mesh:
        #bedinfo and surface info
        md.mesh.vertexonbase = project3d(md, 'vector', np.ones(md.mesh.numberofvertices2d, bool), 'type', 'node', 'layer', 1)
        md.mesh.vertexonsurface = project3d(md, 'vector', np.ones(md.mesh.numberofvertices2d, bool), 'type', 'node', 'layer', md.mesh.numberoflayers)
        md.mesh.vertexonboundary = project3d(md, 'vector', md.mesh.vertexonboundary, 'type', 'node')

        #lat long
        md.mesh.lat = project3d(md, 'vector', md.mesh.lat, 'type', 'node')
        md.mesh.long = project3d(md, 'vector', md.mesh.long, 'type', 'node')
        md.mesh.scale_factor = project3d(md, 'vector', md.mesh.scale_factor, 'type', 'node')

        md.geometry.extrude(md)
        md.friction.extrude(md)
        md.inversion.extrude(md)
        md.smb.extrude(md)
        md.initialization.extrude(md)

        md.flowequation.extrude(md)
        md.stressbalance.extrude(md)
        md.thermal.extrude(md)
        md.masstransport.extrude(md)
        md.levelset.extrude(md)
        md.calving.extrude(md)
        md.frontalforcings.extrude(md)
        md.hydrology.extrude(md)
        md.debris.extrude(md)
        md.solidearth.extrude(md)
        md.dsl.extrude(md)
        md.stochasticforcing.extrude(md)

        #connectivity
        md.mesh.elementconnectivity = np.tile(md.mesh.elementconnectivity, (numlayers - 1, 1))
        md.mesh.elementconnectivity[np.nonzero(md.mesh.elementconnectivity == 0)] = -sys.maxsize - 1
        if not np.isnan(md.mesh.elementconnectivity).all():
            for i in range(1, numlayers - 1):
                connect1 = i * md.mesh.numberofelements2d
                connect2 = (i + 1) * md.mesh.numberofelements2d
                md.mesh.elementconnectivity[connect1:connect2, :] = md.mesh.elementconnectivity[connect1:connect2, :] + md.mesh.numberofelements2d
                md.mesh.elementconnectivity[np.nonzero(md.mesh.elementconnectivity < 0)] = 0

        md.materials.extrude(md)
        md.damage.extrude(md)
        md.mask.extrude(md)
        md.qmu.extrude(md)
        md.basalforcings.extrude(md)
        md.outputdefinition.extrude(md)

        #increase connectivity if less than 25:
        if md.mesh.average_vertex_connectivity <= 25:
            md.mesh.average_vertex_connectivity = 100

        return md
    # }}}

    def collapse(md):  #{{{
        """COLLAPSE - collapses a 3d mesh into a 2d mesh

        This routine collapses a 3d model into a 2d model and collapses all
        the fields of the 3d model by taking their depth-averaged values

        Usage:
            md = collapse(md)

        See also: EXTRUDE, MODELEXTRACT
        """

        # Check that the model is really a 3d model
        if md.mesh.elementtype() != 'Penta':
            raise Exception('collapse error message: only a 3d mesh can be collapsed')

        # Start with changing all the fields from the 3d mesh

        # Dealing with the friction law
        # Drag is limited to nodes that are on the bedrock.
        if md.friction.__class__.__name__ == 'friction':
            md.friction.coefficient = project2d(md, md.friction.coefficient, 1)
            md.friction.p = project2d(md, md.friction.p, 1)
            md.friction.q = project2d(md, md.friction.q, 1)
        elif md.friction.__class__.__name__ == 'frictioncoulomb':
            md.friction.coefficient = project2d(md, md.friction.coefficient, 1)
            md.friction.coefficientcoulomb = project2d(md, md.friction.coefficientcoulomb, 1)
            md.friction.p = project2d(md, md.friction.p, 1)
            md.friction.q = project2d(md, md.friction.q, 1)
        elif md.friction.__class__.__name__ == 'frictionhydro':
            md.friction.q = project2d(md, md.friction.q, 1)
            md.friction.C = project2d(md, md.friction.C, 1)
            md.friction.As = project2d(md, md.friction.As, 1)
            md.friction.effective_pressure = project2d(md, md.friction.effective_pressure, 1)
        elif md.friction.__class__.__name__ == 'frictionwaterlayer':
            md.friction.coefficient = project2d(md, md.friction.coefficient, 1)
            md.friction.p = project2d(md, md.friction.p, 1)
            md.friction.q = project2d(md, md.friction.q, 1)
            md.friction.water_layer = project2d(md, md.friction.water_layer, 1)
        elif md.friction.__class__.__name__ == 'frictionweertman':
            md.friction.C = project2d(md, md.friction.C, 1)
            md.friction.m = project2d(md, md.friction.m, 1)
        elif md.friction.__class__.__name__ == 'frictionweertmantemp':
            md.friction.C = project2d(md, md.friction.C, 1)
            md.friction.m = project2d(md, md.friction.m, 1)
        else:
            print('friction type not supported')

        # Observations
        if not np.isnan(md.inversion.vx_obs).all():
            md.inversion.vx_obs = project2d(md, md.inversion.vx_obs, md.mesh.numberoflayers)
        if not np.isnan(md.inversion.vy_obs).all():
            md.inversion.vy_obs = project2d(md, md.inversion.vy_obs, md.mesh.numberoflayers)
        if not np.isnan(md.inversion.vel_obs).all():
            md.inversion.vel_obs = project2d(md, md.inversion.vel_obs, md.mesh.numberoflayers)
        if not np.isnan(md.inversion.thickness_obs).all():
            md.inversion.thickness_obs = project2d(md, md.inversion.thickness_obs, md.mesh.numberoflayers)
        if not np.isnan(md.inversion.cost_functions_coefficients).all():
            md.inversion.cost_functions_coefficients = project2d(md, md.inversion.cost_functions_coefficients, md.mesh.numberoflayers)
        if isinstance(md.inversion.min_parameters, np.ndarray) and md.inversion.min_parameters.size > 1:
            md.inversion.min_parameters = project2d(md, md.inversion.min_parameters, md.mesh.numberoflayers)
        if isinstance(md.inversion.max_parameters, np.ndarray) and md.inversion.max_parameters.size > 1:
            md.inversion.max_parameters = project2d(md, md.inversion.max_parameters, md.mesh.numberoflayers)
        if md.smb.__class__.__name__ == 'SMBforcing' and not np.isnan(md.smb.mass_balance).all():
            md.smb.mass_balance = project2d(md, md.smb.mass_balance, md.mesh.numberoflayers)
        elif md.smb.__class__.__name__ == 'SMBhenning' and not np.isnan(md.smb.smbref).all():
            md.smb.smbref = project2d(md, md.smb.smbref, md.mesh.numberoflayers)

        # Results
        if not np.isnan(md.initialization.vx).all():
            md.initialization.vx = DepthAverage(md, md.initialization.vx)
        if not np.isnan(md.initialization.vy).all():
            md.initialization.vy = DepthAverage(md, md.initialization.vy)
        if not np.isnan(md.initialization.vz).all():
            md.initialization.vz = DepthAverage(md, md.initialization.vz)
        if not np.isnan(md.initialization.vel).all():
            md.initialization.vel = DepthAverage(md, md.initialization.vel)
        if not np.isnan(md.initialization.temperature).all():
            md.initialization.temperature = DepthAverage(md, md.initialization.temperature)
        if not np.isnan(md.initialization.pressure).all():
            md.initialization.pressure = project2d(md, md.initialization.pressure, 1)
        if not np.isnan(md.initialization.sediment_head).all():
            md.initialization.sediment_head = project2d(md, md.initialization.sediment_head, 1)
        if not np.isnan(md.initialization.epl_head).all():
            md.initialization.epl_head = project2d(md, md.initialization.epl_head, 1)
        if not np.isnan(md.initialization.epl_thickness).all():
            md.initialization.epl_thickness = project2d(md, md.initialization.epl_thickness, 1)
        if not np.isnan(md.initialization.waterfraction).all():
            md.initialization.waterfraction = project2d(md, md.initialization.waterfraction, 1)
        if not np.isnan(md.initialization.watercolumn).all():
            md.initialization.watercolumn = project2d(md, md.initialization.watercolumn, 1)
        if not np.isnan(md.initialization.debris).all():
            md.initialization.debris = project2d(md, md.initialization.debris, 1)

        # elementstype
        if not np.isnan(md.flowequation.element_equation).all():
            md.flowequation.element_equation = project2d(md, md.flowequation.element_equation, 1)
            md.flowequation.vertex_equation = project2d(md, md.flowequation.vertex_equation, 1)
            md.flowequation.borderSSA = project2d(md, md.flowequation.borderSSA, 1)
            md.flowequation.borderHO = project2d(md, md.flowequation.borderHO, 1)
            md.flowequation.borderFS = project2d(md, md.flowequation.borderFS, 1)

        # Boundary conditions
        md.stressbalance.spcvx = project2d(md, md.stressbalance.spcvx, md.mesh.numberoflayers)
        md.stressbalance.spcvy = project2d(md, md.stressbalance.spcvy, md.mesh.numberoflayers)
        md.stressbalance.spcvz = project2d(md, md.stressbalance.spcvz, md.mesh.numberoflayers)
        md.stressbalance.referential = project2d(md, md.stressbalance.referential, md.mesh.numberoflayers)
        md.stressbalance.loadingforce = project2d(md, md.stressbalance.loadingforce, md.mesh.numberoflayers)
        # TODO:
        # - Check if md.mesh.numberoflayershould really be offset by 1.
        # - Find out why md.masstransport.spcthickness is not offset, but the
        #   other fields are.
        # - If offset is required, figure out if it can be abstarcted away to
        #   another part of the API.
        if np.size(md.masstransport.spcthickness) > 1:
            md.masstransport.spcthickness = project2d(md, md.masstransport.spcthickness, md.mesh.numberoflayers)
        if np.size(md.damage.spcdamage) > 1:  # and not np.isnan(md.damage.spcdamage).all():
            md.damage.spcdamage = project2d(md, md.damage.spcdamage, md.mesh.numberoflayers - 1)
        if np.size(md.levelset.spclevelset) > 1:
            md.levelset.spclevelset = project2d(md, md.levelset.spclevelset, md.mesh.numberoflayers - 1)
        md.thermal.spctemperature = project2d(md, md.thermal.spctemperature, md.mesh.numberoflayers - 1)

        # Hydrologydc variables
        if md.hydrology.__class__.__name__ == 'hydrologydc':
            # md.hydrology.spcsediment_head = project2d(md, md.hydrology.spcsediment_head, 1)
            # md.hydrology.mask_eplactive_node = project2d(md, md.hydrology.mask_eplactive_node, 1)
            # md.hydrology.sediment_transmitivity = project2d(md, md.hydrology.sediment_transmitivity, 1)
            # md.hydrology.basal_moulin_input = project2d(md, md.hydrology.basal_moulin_input, 1)
            # if md.hydrology.isefficientlayer == 1:
            #     md.hydrology.spcepl_head = project2d(md, md.hydrology.spcepl_head, 1)
            hydrofields = md.hydrology.__dict__.keys()
            for field in hydrofields:
                try:
                    isvector = np.logical_or(np.shape(md.hydrology.__dict__[field])[0] == md.mesh.numberofelements,
                                             np.shape(md.hydrology.__dict__[field])[0] == md.mesh.numberofvertices)
                except IndexError:
                    isvector = False
                #we collapse only fields that are vertices or element based
                if isvector:
                    md.hydrology.__dict__[field] = project2d(md, md.hydrology.__dict__[field], 1)

        # Materials
        md.materials.rheology_B = DepthAverage(md, md.materials.rheology_B)
        md.materials.rheology_n = project2d(md, md.materials.rheology_n, 1)

        # dsl
        if np.size(md.dsl.sea_surface_height_above_geoid) > 1:
            md.dsl.sea_surface_height_above_geoid = project2d(md, md.dsl.sea_surface_height_above_geoid, 1)
        if np.size(md.dsl.sea_water_pressure_at_sea_floor) > 1:
            md.dsl.sea_water_pressure_at_sea_floor = project2d(md, md.dsl.sea_water_pressure_at_sea_floor, 1)

        # Damage
        if md.damage.isdamage:
            md.damage.D = DepthAverage(md, md.damage.D)

        # Special for thermal modeling
        if not np.isnan(md.basalforcings.groundedice_melting_rate).all():
            md.basalforcings.groundedice_melting_rate = project2d(md, md.basalforcings.groundedice_melting_rate, 1)
        if hasattr(md.basalforcings, 'floatingice_melting_rate') and not np.isnan(md.basalforcings.floatingice_melting_rate).all():
            md.basalforcings.floatingice_melting_rate = project2d(md, md.basalforcings.floatingice_melting_rate, 1)
        md.basalforcings.geothermalflux = project2d(md, md.basalforcings.geothermalflux, 1)  #bedrock only gets geothermal flux

        if hasattr(md.calving, 'coeff') and not np.isnan(md.calving.coeff).all():
            md.calving.coeff = project2d(md, md.calving.coeff, 1)
        if hasattr(md.frontalforcings, 'meltingrate') and not np.isnan(md.frontalforcings.meltingrate).all():
            md.frontalforcings.meltingrate = project2d(md, md.frontalforcings.meltingrate, 1)

        # Update of connectivity matrix
        md.mesh.average_vertex_connectivity = 25

        # Collapse the mesh
        nodes2d = md.mesh.numberofvertices2d
        elements2d = md.mesh.numberofelements2d

        # Parameters
        md.geometry.surface = project2d(md, md.geometry.surface, 1)
        md.geometry.thickness = project2d(md, md.geometry.thickness, 1)
        md.geometry.base = project2d(md, md.geometry.base, 1)
        if not np.isnan(md.geometry.bed).all():
            md.geometry.bed = project2d(md, md.geometry.bed, 1)
        if not np.isnan(md.mask.ocean_levelset).all():
            md.mask.ocean_levelset = project2d(md, md.mask.ocean_levelset, 1)
        if not np.isnan(md.mask.ice_levelset).all():
            md.mask.ice_levelset = project2d(md, md.mask.ice_levelset, 1)

        # Lat/long
        if np.size(md.mesh.lat) == md.mesh.numberofvertices:
            md.mesh.lat = project2d(md, md.mesh.lat, 1)
        if np.size(md.mesh.long) == md.mesh.numberofvertices:
            md.mesh.long = project2d(md, md.mesh.long, 1)

        # OutputDefinitions
        if md.outputdefinition.definitions:
            for solutionfield, field in list(md.outputdefinition.__dict__.items()):
                if isinstance(field, list):
                    # Get each definition
                    for i, fieldi in enumerate(field):
                        if fieldi:
                            fieldr = getattr(md.outputdefinition, solutionfield)[i]
                            # Get subfields
                            for solutionsubfield, subfield in list(fieldi.__dict__.items()):
                                if np.size(subfield) == md.mesh.numberofvertices:
                                    setattr(fieldr, solutionsubfield, project2d(md, subfield, 1))
                                elif np.size(subfield) == md.mesh.numberofelements:
                                    setattr(fieldr, solutionsubfield, project2d(md, subfield, 1))

        # Initialize 2d mesh
        mesh = mesh2d()
        mesh.x = md.mesh.x2d
        mesh.y = md.mesh.y2d
        mesh.numberofvertices = md.mesh.numberofvertices2d
        mesh.numberofelements = md.mesh.numberofelements2d
        mesh.elements = md.mesh.elements2d
        # if not np.isnan(md.mesh.vertexonboundary).all():
        #     mesh.vertexonboundary = project2d(md, md.mesh.vertexonboundary, 1)
        # if not np.isnan(md.mesh.elementconnectivity).all():
        #     mesh.elementconnectivity = project2d(md, md.mesh.elementconnectivity, 1)
        if np.size(md.mesh.lat) == md.mesh.numberofvertices:
            mesh.lat = project2d(md, md.mesh.lat, 1)
        if np.size(md.mesh.long) == md.mesh.numberofvertices:
            mesh.long = project2d(md, md.mesh.long, 1)
        mesh.epsg = md.mesh.epsg
        if np.size(md.mesh.scale_factor) == md.mesh.numberofvertices:
            mesh.scale_factor = project2d(md, md.mesh.scale_factor, 1)
        if hasattr(md.mesh, 'vertexonboundary') and not np.isnan(md.mesh.vertexonboundary).all():
            mesh.vertexonboundary = project2d(md, md.mesh.vertexonboundary, 1)
        if hasattr(md.mesh, 'elementonboundary') and not np.isnan(md.mesh.elementonboundary).all():
            mesh.elementonboundary = project2d(md, md.mesh.elementonboundary, 1)
        md.mesh = mesh
        md.mesh.vertexconnectivity = NodeConnectivity(md.mesh.elements, md.mesh.numberofvertices)
        md.mesh.elementconnectivity = ElementConnectivity(md.mesh.elements, md.mesh.vertexconnectivity)
        md.mesh.segments = contourenvelope(md.mesh)

        return md
    # }}}
