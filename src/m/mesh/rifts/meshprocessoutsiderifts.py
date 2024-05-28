import numpy as np

from ContourToMesh import ContourToMesh
from ElementsFromEdge import ElementsFromEdge
import MatlabFuncs as m


def meshprocessoutsiderifts(md, domainoutline):
    """MESHPROCESSOUTSIDERIFTS - process rifts when they touch the domain outline

    Usage:
        md = meshprocessoutsiderifts(md, domain)
    """

    #go through rifts, and figure out which ones touch the domain outline
    for rift in md.rifts.riftstruct:

        #first, flag nodes that belong to the domain outline
        flags = ContourToMesh(md.mesh.elements, md.mesh.x, md.mesh.y, domainoutline, 'node', 0)

        tips = rift.tips
        outsidetips = tips[np.nonzero(flags[rift.tips - 1])[0]]

        #we have found outsidetips, tips that touch the domain outline. go through them
        for tip in outsidetips:
            #find tip in the segments, take first segment (there should be 2) that holds tip,
            #and node_connected_to_tip is the other node on this segment:
            tipindex = np.nonzero(rift.segments[:, 0] == tip)[0]
            if tipindex:
                tipindex = tipindex[0]
                node_connected_to_tip = rift.segments[tipindex, 1]
            else:
                tipindex = np.nonzero(rift.segments[:, 1] == tip)[0]
                tipindex = tipindex[0]
                node_connected_to_tip = rift.segments[tipindex, 1]

            #ok, we have the tip node, and the first node connected to it, on the rift. Now,
            #identify all the elements that are connected to the tip, and that are on the same
            #side of the rift.
            A = tip
            B = node_connected_to_tip

            elements = np.empty(0, int)

            while flags(B):  #as long as B does not belong to the domain outline, keep looking.
                #detect elements on edge A, B:
                edgeelements = ElementsFromEdge(md.mesh.elements, A, B)
                #rule out those we already detected
                already_detected = m.ismember(edgeelements, elements)
                nextelement = edgeelements(np.nonzero(np.logical_not(already_detected))[0])
                #add new detected element to the list of elements we are looking for.
                elements = np.concatenate((elements, nextelement))
                #new B:
                B = md.mesh.elements[nextelement - 1, np.nonzero(np.logical_not(m.ismember(md.mesh.elements[nextelement - 1, :], np.array([A, B]))))[0]]

            #take the list of elements on one side of the rift that connect to the tip,
            #and duplicate the tip on them, so as to open the rift to the outside.
            num = np.size(md.mesh.x) + 1
            md.mesh.x = np.concatenate((md.mesh.x, md.mesh.x[tip]))
            md.mesh.y = np.concatenate((md.mesh.y, md.mesh.y[tip]))
            md.mesh.numberofvertices = num

            #replace tip in elements
            newelements = md.mesh.elements[elements - 1, :]
            pos = np.nonzero(newelements == tip)[0]
            newelements[pos] = num
            md.mesh.elements[elements - 1, :] = newelements
            rift.tips = np.concatenate((rift.tips, num))

            #deal with segments
            tipsegments = np.nonzero(np.logical_or(md.mesh.segments[:, 0] == tip, md.mesh.segments[:, 1] == tip))[0]
            for segment_index in tipsegments:
                pos = np.nonzero(md.mesh.segments[segment_index, 0:2] != tip)[0]
                other_node = md.mesh.segments[segment_index, pos]
                if not isconnected(md.mesh.elements, other_node, tip):
                    pos = np.nonzero(md.mesh.segments[segment_index, 0:2] == tip)[0]
                    md.mesh.segments[segment_index, pos] = num

    #Fill in rest of fields:
    md.mesh.numberofelements = np.size(md.mesh.elements, axis=0)
    md.mesh.numberofvertices = np.size(md.mesh.x)
    md.mesh.vertexonboundary = np.zeros(np.size(md.mesh.x), int)
    md.mesh.vertexonboundary[md.mesh.segments[:, 0:2] - 1] = 1
    md.rifts.numrifts = np.length(md.rifts.riftstruct)

    return md


def isconnected(elements, A, B):  #{{{
    """ISCONNECTED: are two nodes connected by a triangulation?

    Usage:
        flag = isconnected(elements, A, B)
    """

    elements = ElementsFromEdge(elements, A, B)
    if not elements:
        flag = 0
    else:
        flag = 1

    return flag
    # }}}
