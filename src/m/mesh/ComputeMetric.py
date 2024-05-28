import numpy as np


def ComputeMetric(hessian, scale, epsilon, hmin, hmax, pos):
    """
    COMPUTEMETRIC - compute metric from an Hessian

       Usage:
          metric = ComputeMetric(hessian, scale, epsilon, hmin, hmax, pos)
          pos is contains the positions where the metric is wished to be maximized (water?)

       Example:
          metric = ComputeMetric(hessian, 2 / 9, 1.0e-1, 100, 1.0e5, [])
    """

    #first, find the eigen values of each line of H = [hessian(i, 1) hessian(i, 2); hessian(i, 2) hessian(i, 3)]
    a = hessian[:, 0]
    b = hessian[:, 1]
    d = hessian[:, 2]
    lambda1 = 0.5 * ((a + d) + np.sqrt(4. * b**2 + (a - d)**2))
    lambda2 = 0.5 * ((a + d) - np.sqrt(4. * b**2 + (a - d)**2))
    pos1 = np.nonzero(lambda1 == 0.)[0]
    pos2 = np.nonzero(lambda2 == 0.)[0]
    pos3 = np.nonzero(np.logical_and(b == 0., lambda1 == lambda2))[0]

    #Modify the eigen values to control the shape of the elements
    lambda1 = np.minimum(np.maximum(np.abs(lambda1) * scale / epsilon, 1. / hmax**2), 1. / hmin**2)
    lambda2 = np.minimum(np.maximum(np.abs(lambda2) * scale / epsilon, 1. / hmax**2), 1. / hmin**2)

    #compute eigen vectors
    norm1 = np.sqrt(8. * b**2 + 2. * (d - a)**2 + 2. * (d - a) * np.sqrt((a - d)**2 + 4. * b**2))
    v1x = 2. * b / norm1
    v1y = ((d - a) + np.sqrt((a - d)**2 + 4. * b**2)) / norm1
    norm2 = np.sqrt(8. * b**2 + 2. * (d - a)**2 - 2. * (d - a) * np.sqrt((a - d)**2 + 4. * b**2))
    v2x = 2. * b / norm2
    v2y = ((d - a) - np.sqrt((a - d)**2 + 4. * b**2)) / norm2

    v1x[pos3] = 1.
    v1y[pos3] = 0.
    v2x[pos3] = 0.
    v2y[pos3] = 1.

    #Compute new metric (for each node M = V * Lambda * V^-1)

    metric = np.vstack((((v1x * v2y - v1y * v2x)**(-1) * (lambda1 * v2y * v1x - lambda2 * v1y * v2x)).reshape(-1, ),
                        ((v1x * v2y - v1y * v2x)**(-1) * (lambda1 * v1y * v2y - lambda2 * v1y * v2y)).reshape(-1, ),
                        ((v1x * v2y - v1y * v2x)**(-1) * (-lambda1 * v2x * v1y + lambda2 * v1x * v2y)).reshape(-1, ))).T

    #some corrections for 0 eigen values
    metric[pos1, :] = np.tile(np.array([[1. / hmax**2, 0., 1. / hmax**2]]), (np.size(pos1), 1))
    metric[pos2, :] = np.tile(np.array([[1. / hmax**2, 0., 1. / hmax**2]]), (np.size(pos2), 1))

    #take care of water elements
    metric[pos, :] = np.tile(np.array([[1. / hmax**2, 0., 1. / hmax**2]]), (np.size(pos), 1))

    #take care of NaNs if any (use Numpy eig in a loop)
    pos = np.nonzero(np.isnan(metric))[0]
    if np.size(pos):
        print((" %i NaN found in the metric. Use Numpy routine..." % np.size(pos)))
        for posi in pos:
            H = np.array([[hessian[posi, 0], hessian[posi, 1]], [hessian[posi, 1], hessian[posi, 2]]])
            [v, u] = np.linalg.eig(H)
            v = np.diag(v)
            lambda1 = v[0, 0]
            lambda2 = v[1, 1]
            v[0, 0] = np.minimum(np.maximum(np.abs(lambda1) * scale / epsilon, 1. / hmax**2), 1. / hmin**2)
            v[1, 1] = np.minimum(np.maximum(np.abs(lambda2) * scale / epsilon, 1. / hmax**2), 1. / hmin**2)

            metricTria = np.dot(np.dot(u, v), np.linalg.inv(u))
            metric[posi, :] = np.array([metricTria[0, 0], metricTria[0, 1], metricTria[1, 1]])

    if np.any(np.isnan(metric)):
        raise RunTimeError("ComputeMetric error message: NaN in the metric despite our efforts...")

    return metric
