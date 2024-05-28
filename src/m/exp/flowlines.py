import numpy as np
from pairoptions import pairoptions
try:
    import matplotlib.tri as tri
except ImportError:
    print("could not import matplotlib which is needed for triangulation")


def flowlines(index, x, y, u, v, x0, y0, *args):
    #FLOWLINES - compute flowlines from a given set of seed points
    #
    #   Usage:
    #      flowpath = flowlines(index, x, y, u, v, x0, y0, varargin)
    #
    #   the velocity field is given by the couple (u, v) and the coordinates
    #   of the seed points are (x0, y0). One can use one or several seed
    #   points
    #
    #   Example:
    #      flowpath = flowlines(md.mesh.elements, md.mesh.x, md.mesh.y, md.initialization.vx, md.initialization.vy, x0, y0)
    #
    #   Options:
    # - 'maxiter':   how many steps upstream and downstream of the seed points (default: 200)
    # - 'precision': division of each segment (higer precision increases number of segments, default: 1)
    # - 'downstream':flow line upstream of the seed points (default: 1)
    # - 'upstream':  flow line upstream of the seed points (default: 1)

    #check input
    if (not len(x) == len(y) == len(u) == len(v)):
        raise IOError('flowlines error message: x, y, u and v must have the same length')

    if len(x) < 3:
        raise IOError('flowlines error message: at least one element is required')

    if len(x0) != len(y0):
        raise IOError('flowlines error message: x0 and y0 do not have the same length')

    #check if we have matlab indexing and translate
    if np.nanmin(index) > 0:
        index -= 1

    #process options
    options = pairoptions(*args)
    maxiter = options.getfieldvalue('maxiter', 200)
    precision = options.getfieldvalue('precision', 1)

    #Squeeze the result to fix shape issues
    u = np.squeeze(u)
    v = np.squeeze(v)

    #Create triangulation once for all and check seed points
    trep = tri.Triangulation(x, y, index)
    trifinder = trep.get_trifinder()
    tria = trifinder(x0, y0)
    pos = np.where(tria < 0)[0]
    if len(pos) > 0:
        x0 = np.delete(x0, pos)
        y0 = np.delete(y0, pos)

    #initialize other variables
    N = len(x0)
    flowpath = {}
    flowpath['x'] = list.copy(x0)
    flowpath['y'] = list.copy(y0)
    flowpath['name'] = 'flowline{}'.format(0)
    for i in np.arange(1, len(flowpath['x'])):
        flowpath['name'] = np.append(flowpath['name'], 'flowline{}'.format(i))

    #get avegared length of each element
    dist0_1 = np.sqrt((x[index[:, 0]] - x[index[:, 1]])**2 + (y[index[:, 0]] - y[index[:, 1]])**2)
    dist0_2 = np.sqrt((x[index[:, 0]] - x[index[:, 2]])**2 + (y[index[:, 0]] - y[index[:, 2]])**2)
    dist1_2 = np.sqrt((x[index[:, 1]] - x[index[:, 2]])**2 + (y[index[:, 1]] - y[index[:, 2]])**2)
    length_tria = 1 / 3 * (dist0_1 + dist0_2 + dist1_2)

    #take velocity for each element
    u = np.nanmean(u[index], axis=1)
    v = np.nanmean(v[index], axis=1)

    for flowdirection in ['downstream', 'upstream']:
        print('Dealing with the {} flowlines'.format(flowdirection))
    #initialization:
        counter = 1
        treatdirection = options.getfieldvalue(flowdirection, 1)
        done = np.zeros(N)
        queue = []
        X = np.array(x0)
        Y = np.array(y0)

        if treatdirection:
            if flowdirection == 'upstream':
                flowindex = 0
            elif flowdirection == 'downstream':
                flowindex = -1

            while not all(done):
                #find current triangle
                queue = np.where(done == 0)[0]
                tria = trifinder(X[queue], Y[queue])
                #check that the point are actually inside a triangle of the mesh
                outsiders = np.where(tria < 0)[0]
                print('outsider is {} with queue {}'.format(outsiders, queue))
                if np.size(outsiders) > 0:
                    for outsider in outsiders:
                        flowpath['x'][queue[outsider]] = np.delete(flowpath['x'][queue[outsider]][:], flowindex)
                        flowpath['y'][queue[outsider]] = np.delete(flowpath['y'][queue[outsider]][:], flowindex)

                    done[queue[outsiders]] = 1

                    tria = np.delete(tria, outsiders)
                    queue = np.delete(queue, outsiders)

                #velocity of the current triangle and norm it
                if flowdirection == 'upstream':
                    ut = -u[tria]
                    vt = -v[tria]
                if flowdirection == 'downstream':
                    ut = u[tria]
                    vt = v[tria]
                normv = np.sqrt(ut**2 + vt**2)
                normv[np.where(normv < 1.0e-10)] = 1.0e-10
                ut = ut / normv
                vt = vt / normv

                if counter > maxiter:
                    print('Maximum number of iterations ({}) reached while going {}'.format(maxiter, flowdirection))
                    break

                counter += 1

                #remove stagnant point
                stagnants = np.where(np.logical_and(ut == 0, vt == 0))
                done[queue[stagnants]] = 1
                #build next point
                for i in np.arange(0, len(queue)):
                    if np.size(flowpath['x'][queue[i]]) == 1:
                        X[queue[i]] = flowpath['x'][queue[i]] + ut[i] * length_tria[tria[i]] / precision
                        Y[queue[i]] = flowpath['y'][queue[i]] + vt[i] * length_tria[tria[i]] / precision
                    else:
                        X[queue[i]] = flowpath['x'][queue[i]][flowindex] + ut[i] * length_tria[tria[i]] / precision
                        Y[queue[i]] = flowpath['y'][queue[i]][flowindex] + vt[i] * length_tria[tria[i]] / precision
                    if flowdirection == 'upstream':
                        flowpath['x'][queue[i]] = np.append(X[queue[i]], flowpath['x'][queue[i]])
                        flowpath['y'][queue[i]] = np.append(Y[queue[i]], flowpath['y'][queue[i]])
                    elif flowdirection == 'downstream':
                        flowpath['x'][queue[i]] = np.append(flowpath['x'][queue[i]], X[queue[i]])
                        flowpath['y'][queue[i]] = np.append(flowpath['y'][queue[i]], Y[queue[i]])

    return flowpath
