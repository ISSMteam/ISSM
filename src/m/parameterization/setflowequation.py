import numpy as np
from model import model
from pairoptions import pairoptions
from FlagElements import FlagElements


def setflowequation(md, *args):
    """SETFLOWEQUATION - associate a solution type to each element

    This routine works like plotmodel: it works with an even number of inputs
    'SIA', 'SSA', 'HO', 'L1L2', 'MOLHO', 'FS' and 'fill' are the possible 
    options that must be followed by the corresponding exp file or flags list. 
    It can either be a domain file (argus type, .exp extension), or an array of 
    element flags.
    If user wants every element outside the domain to be setflowequationd, add 
    '~' to the name of the domain file (ex: '~HO.exp') an empty string '' will 
    be considered as an empty domain a string 'all' will be considered as the 
    entire domain.
    You can specify the type of coupling, 'penalties' or 'tiling', to use with 
    the input 'coupling'.

    Usage:
        md = setflowequation(md, varargin)

    Example:
        md = setflowequation(md, 'HO', 'HO.exp', fill', 'SIA', 'coupling', 'tiling')
    """

    #some checks on list of arguments
    if not isinstance(md, model) or not len(args):
        raise TypeError("setflowequation error message")

    #process options
    options = pairoptions(*args)
    #    options = deleteduplicates(options, 1)

    #Find_out what kind of coupling to use
    coupling_method = options.getfieldvalue('coupling', 'tiling')
    if coupling_method not in ['tiling', 'penalties']:
        raise TypeError("coupling type can only be: tiling or penalties")

    #recover elements distribution
    SIAflag = FlagElements(md, options.getfieldvalue('SIA', ''))
    SSAflag = FlagElements(md, options.getfieldvalue('SSA', ''))
    HOflag = FlagElements(md, options.getfieldvalue('HO', ''))
    L1L2flag = FlagElements(md, options.getfieldvalue('L1L2', ''))
    MOLHOflag = FlagElements(md, options.getfieldvalue('MOLHO', ''))
    FSflag = FlagElements(md, options.getfieldvalue('FS', ''))
    filltype = options.getfieldvalue('fill', 'none')

    #Flag the elements that have not been flagged as filltype
    if 'SIA' in filltype:
        SIAflag = ~SSAflag & ~HOflag
    elif 'SSA' in filltype:
        SSAflag = ~SIAflag & ~HOflag & ~FSflag
    elif 'HO' in filltype:
        HOflag = ~SIAflag & ~SSAflag & ~FSflag
    #check that each element has at least one flag
    if not any(SIAflag + SSAflag + L1L2flag + MOLHOflag + HOflag + FSflag):
        raise TypeError("elements type not assigned, supported models are 'SIA', 'SSA', 'HO' and 'FS'")

    #check that each element has only one flag
    if any(SIAflag + SSAflag + L1L2flag + MOLHOflag + HOflag + FSflag > 1):
        print('Warning: setflowequation.py: some elements have several types, higher order type is used for them')
        SIAflag[np.where(np.logical_and(SIAflag, SSAflag))] = False
        SIAflag[np.where(np.logical_and(SIAflag, HOflag))] = False
        SSAflag[np.where(np.logical_and(SSAflag, HOflag))] = False

        #check that L1L2 and MOLHO is not coupled to any other model for now
        if any(L1L2flag) and any(SIAflag + SSAflag + HOflag + FSflag):
            raise TypeError('L1L2 cannot be coupled to any other model')
        if any(MOLHOflag) and any(SIAflag + SSAflag + HOflag + FSflag):
            raise TypeError('MOLHO cannot be coupled to any other model')

        #Check that no HO or FS for 2d mesh
        if md.mesh.domaintype == '2Dhorizontal':
            if any(FSflag + HOflag):
                raise TypeError('FS and HO elements not allowed in 2d mesh, extrude it first')

    #FS can only be used alone for now:
    if any(FSflag) and any(SIAflag):
        raise TypeError("FS cannot be used with any other model for now, put FS everywhere")

    #Initialize node fields
    nodeonSIA = np.zeros(md.mesh.numberofvertices, bool)
    nodeonSIA[md.mesh.elements[np.where(SIAflag), :] - 1] = True
    nodeonSSA = np.zeros(md.mesh.numberofvertices, bool)
    nodeonSSA[md.mesh.elements[np.where(SSAflag), :] - 1] = True
    nodeonL1L2 = np.zeros(md.mesh.numberofvertices, bool)
    nodeonL1L2[md.mesh.elements[np.where(L1L2flag), :] - 1] = True
    nodeonMOLHO = np.zeros(md.mesh.numberofvertices, bool)
    nodeonMOLHO[md.mesh.elements[np.where(MOLHOflag), :] - 1] = True
    nodeonHO = np.zeros(md.mesh.numberofvertices, bool)
    nodeonHO[md.mesh.elements[np.where(HOflag), :] - 1] = True
    nodeonFS = np.zeros(md.mesh.numberofvertices, bool)
    noneflag = np.zeros(md.mesh.numberofelements, bool)

    #First modify FSflag to get rid of elements contrained everywhere (spc + border with HO or SSA)
    if any(FSflag):
        fullspcnodes = np.logical_or(~np.isnan(md.stressbalance.spcvx) & ~np.isnan(md.stressbalance.spcvy) & ~np.isnan(md.stressbalance.spcvz), np.logical_and(nodeonHO, nodeonFS))  #find all the nodes on the boundary of the domain without icefront
        fullspcelems = np.sum(fullspcnodes[md.mesh.elements - 1], axis=1) == 6  #find all the nodes on the boundary of the domain without icefront
        FSflag[np.where(fullspcelems.reshape(-1))] = False
        nodeonFS[md.mesh.elements[np.where(FSflag), :] - 1] = True

    #Then complete with NoneApproximation or the other model used if there is no FS
    if any(FSflag):
        if any(HOflag):  #fill with HO
            HOflag[~FSflag] = True
            nodeonHO[md.mesh.elements[np.where(HOflag), :] - 1] = True
        elif any(SSAflag):  #fill with SSA
            SSAflag[~FSflag] = True
            nodeonSSA[md.mesh.elements[np.where(SSAflag), :] - 1] = True
        else:  #fill with none
            noneflag[np.where(~FSflag)] = True

    #Now take care of the coupling between SSA and HO
    if coupling_method not in ['penalties']:
        md.stressbalance.vertex_pairing = np.array([])
    nodeonSSAHO = np.zeros(md.mesh.numberofvertices, bool)
    nodeonHOFS = np.zeros(md.mesh.numberofvertices, bool)
    nodeonSSAFS = np.zeros(md.mesh.numberofvertices, bool)
    SSAHOflag = np.zeros(md.mesh.numberofelements, bool)
    SSAFSflag = np.zeros(md.mesh.numberofelements, bool)
    HOFSflag = np.zeros(md.mesh.numberofelements, bool)
    if coupling_method == 'penalties':
        #Create the border nodes between HO and SSA and extrude them
        numnodes2d = md.mesh.numberofvertices2d
        numlayers = md.mesh.numberoflayers
        bordernodes2d = np.where(np.logical_and(nodeonHO[0:numnodes2d], nodeonSSA[0:numnodes2d]))[0] + 1  #Nodes connected to two different types of elements

    #initialize and fill in penalties structure
        if np.all(np.logical_not(np.isnan(bordernodes2d))):
            penalties = np.zeros((0, 2))
            for i in range(1, numlayers):
                penalties = np.vstack((penalties, np.vstack((bordernodes2d, bordernodes2d + md.mesh.numberofvertices2d * (i))).T))
            md.stressbalance.vertex_pairing = penalties

    elif coupling_method == 'tiling':
        if any(SSAflag) and any(HOflag):  #coupling SSA HO
            #Find node at the border
            nodeonSSAHO[np.where(np.logical_and(nodeonSSA, nodeonHO))] = True
            #SSA elements in contact with this layer become SSAHO elements
            matrixelements = nodeonSSAHO[md.mesh.elements - 1]
            commonelements = np.sum(matrixelements, axis=1) != 0
            commonelements[np.where(HOflag)] = False  #only one layer: the elements previously in SSA
            SSAflag[np.where(commonelements)] = False  #these elements are now SSAHOelements
            SSAHOflag[np.where(commonelements)] = True
            nodeonSSA[:] = False
            nodeonSSA[md.mesh.elements[np.where(SSAflag), :] - 1] = True

            #rule out elements that don't touch the 2 boundaries
            pos = np.where(SSAHOflag)[0]
            elist = np.zeros(np.size(pos), dtype=int)
            elist = elist + np.sum(nodeonSSA[md.mesh.elements[pos, :] - 1], axis=1).astype(bool)
            elist = elist - np.sum(nodeonHO[md.mesh.elements[pos, :] - 1], axis=1).astype(bool)
            pos1 = np.where(elist == 1)[0]
            SSAflag[pos[pos1]] = True
            SSAHOflag[pos[pos1]] = False
            pos2 = np.where(elist == -1)[0]
            HOflag[pos[pos2]] = True
            SSAHOflag[pos[pos2]] = False

            #Recompute nodes associated to these elements
            nodeonSSA[:] = False
            nodeonSSA[md.mesh.elements[np.where(SSAflag), :] - 1] = True
            nodeonHO[:] = False
            nodeonHO[md.mesh.elements[np.where(HOflag), :] - 1] = True
            nodeonSSAHO[:] = False
            nodeonSSAHO[md.mesh.elements[np.where(SSAHOflag), :] - 1] = True

        elif any(HOflag) and any(FSflag):  #coupling HO FS
            #Find node at the border
            nodeonHOFS[np.where(np.logical_and(nodeonHO, nodeonFS))] = True
            #FS elements in contact with this layer become HOFS elements
            matrixelements = nodeonHOFS[md.mesh.elements - 1]
            commonelements = np.sum(matrixelements, axis=1) != 0
            commonelements[np.where(HOflag)] = False  #only one layer: the elements previously in SSA
            FSflag[np.where(commonelements)] = False  #these elements are now SSAHOelements
            HOFSflag[np.where(commonelements)] = True
            nodeonFS = np.zeros(md.mesh.numberofvertices, bool)
            nodeonFS[md.mesh.elements[np.where(FSflag), :] - 1] = True

            #rule out elements that don't touch the 2 boundaries
            pos = np.where(HOFSflag)[0]
            elist = np.zeros(np.size(pos), dtype=int)
            elist = elist + np.sum(nodeonFS[md.mesh.elements[pos, :] - 1], axis=1).astype(bool)
            elist = elist - np.sum(nodeonHO[md.mesh.elements[pos, :] - 1], axis=1).astype(bool)
            pos1 = np.where(elist == 1)[0]
            FSflag[pos[pos1]] = True
            HOFSflag[pos[pos1]] = False
            pos2 = np.where(elist == -1)[0]
            HOflag[pos[pos2]] = True
            HOFSflag[pos[pos2]] = False

            #Recompute nodes associated to these elements
            nodeonFS[:] = False
            nodeonFS[md.mesh.elements[np.where(FSflag), :] - 1] = True
            nodeonHO[:] = False
            nodeonHO[md.mesh.elements[np.where(HOflag), :] - 1] = True
            nodeonHOFS[:] = False
            nodeonHOFS[md.mesh.elements[np.where(HOFSflag), :] - 1] = True
        elif any(FSflag) and any(SSAflag):
            #Find node at the border
            nodeonSSAFS[np.where(np.logical_and(nodeonSSA, nodeonFS))] = True
            #FS elements in contact with this layer become SSAFS elements
            matrixelements = nodeonSSAFS[md.mesh.elements - 1]
            commonelements = np.sum(matrixelements, axis=1) != 0
            commonelements[np.where(SSAflag)] = False  #only one layer: the elements previously in SSA
            FSflag[np.where(commonelements)] = False  #these elements are now SSASSAelements
            SSAFSflag[np.where(commonelements)] = True
            nodeonFS = np.zeros(md.mesh.numberofvertices, bool)
            nodeonFS[md.mesh.elements[np.where(FSflag), :] - 1] = True

            #rule out elements that don't touch the 2 boundaries
            pos = np.where(SSAFSflag)[0]
            elist = np.zeros(np.size(pos), dtype=int)
            elist = elist + np.sum(nodeonSSA[md.mesh.elements[pos, :] - 1], axis=1).astype(bool)
            elist = elist - np.sum(nodeonFS[md.mesh.elements[pos, :] - 1], axis=1).astype(bool)
            pos1 = np.where(elist == 1)[0]
            SSAflag[pos[pos1]] = True
            SSAFSflag[pos[pos1]] = False
            pos2 = np.where(elist == -1)[0]
            FSflag[pos[pos2]] = True
            SSAFSflag[pos[pos2]] = False

            #Recompute nodes associated to these elements
            nodeonSSA[:] = False
            nodeonSSA[md.mesh.elements[np.where(SSAflag), :] - 1] = True
            nodeonFS[:] = False
            nodeonFS[md.mesh.elements[np.where(FSflag), :] - 1] = True
            nodeonSSAFS[:] = False
            nodeonSSAFS[md.mesh.elements[np.where(SSAFSflag), :] - 1] = True

        elif any(FSflag) and any(SIAflag):
            raise TypeError("type of coupling not supported yet")

    #Create SSAHOApproximation where needed
    md.flowequation.element_equation = np.zeros(md.mesh.numberofelements, int)
    md.flowequation.element_equation[np.where(noneflag)] = 0
    md.flowequation.element_equation[np.where(SIAflag)] = 1
    md.flowequation.element_equation[np.where(SSAflag)] = 2
    md.flowequation.element_equation[np.where(L1L2flag)] = 3
    md.flowequation.element_equation[np.where(MOLHOflag)] = 4
    md.flowequation.element_equation[np.where(HOflag)] = 5
    md.flowequation.element_equation[np.where(FSflag)] = 6
    md.flowequation.element_equation[np.where(SSAHOflag)] = 7
    md.flowequation.element_equation[np.where(SSAFSflag)] = 8
    md.flowequation.element_equation[np.where(HOFSflag)] = 9

    #border
    md.flowequation.borderHO = nodeonHO
    md.flowequation.borderSSA = nodeonSSA
    md.flowequation.borderFS = nodeonFS

    #Create vertices_type
    md.flowequation.vertex_equation = np.zeros(md.mesh.numberofvertices, int)
    pos = np.where(nodeonSSA)
    md.flowequation.vertex_equation[pos] = 2
    pos = np.where(nodeonL1L2)
    md.flowequation.vertex_equation[pos] = 3
    pos = np.where(nodeonMOLHO)
    md.flowequation.vertex_equation[pos] = 4
    pos = np.where(nodeonHO)
    md.flowequation.vertex_equation[pos] = 5
    pos = np.where(nodeonFS)
    md.flowequation.vertex_equation[pos] = 6
    #DO SIA LAST! Otherwise spcs might not be set up correctly (SIA should have priority)
    pos = np.where(nodeonSIA)
    md.flowequation.vertex_equation[pos] = 1
    if any(FSflag):
        pos = np.where(np.logical_not(nodeonFS))
        if not (any(HOflag) or any(SSAflag)):
            md.flowequation.vertex_equation[pos] = 0
    pos = np.where(nodeonSSAHO)
    md.flowequation.vertex_equation[pos] = 7
    pos = np.where(nodeonHOFS)
    md.flowequation.vertex_equation[pos] = 8
    pos = np.where(nodeonSSAFS)
    md.flowequation.vertex_equation[pos] = 9

    #figure out solution types
    md.flowequation.isSIA = any(md.flowequation.element_equation == 1)
    md.flowequation.isSSA = any(md.flowequation.element_equation == 2)
    md.flowequation.isL1L2= any(md.flowequation.element_equation == 3)
    md.flowequation.isMOLHO= any(md.flowequation.element_equation == 4)
    md.flowequation.isHO = any(md.flowequation.element_equation == 5)
    md.flowequation.isFS = any(md.flowequation.element_equation == 6)

    return md

    #Check that tiling can work:
    if any(md.flowequation.borderSSA) and any(md.flowequation.borderHO) and any(md.flowequation.borderHO + md.flowequation.borderSSA != 1):
        raise TypeError("error coupling domain too irregular")
    if any(md.flowequation.borderSSA) and any(md.flowequation.borderFS) and any(md.flowequation.borderFS + md.flowequation.borderSSA != 1):
        raise TypeError("error coupling domain too irregular")
    if any(md.flowequation.borderFS) and any(md.flowequation.borderHO) and any(md.flowequation.borderHO + md.flowequation.borderFS != 1):
        raise TypeError("error coupling domain too irregular")

    return md
