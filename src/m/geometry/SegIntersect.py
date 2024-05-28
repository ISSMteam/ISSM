import numpy as np


def SegIntersect(seg1, seg2):
    """
    SEGINTERSECT - test of segments intersection

       return 1 if the two segments intersect
       seg1 = [x1 y1; x2 y2]
       seg2 = [x1 y1; x2 y2]

       Usage:
          bval = SegIntersect(seg1, seg2)
    """

    bval = 1

    xA = seg1[0, 0]
    yA = seg1[0, 1]
    xB = seg1[1, 0]
    yB = seg1[1, 1]
    xC = seg2[0, 0]
    yC = seg2[0, 1]
    xD = seg2[1, 0]
    yD = seg2[1, 1]

    O2A = np.array([xA, yA]) - np.array([xD / 2. + xC / 2., yD / 2. + yC / 2.])
    O2B = np.array([xB, yB]) - np.array([xD / 2. + xC / 2., yD / 2. + yC / 2.])
    O1C = np.array([xC, yC]) - np.array([xA / 2. + xB / 2., yB / 2. + yA / 2.])
    O1D = np.array([xD, yD]) - np.array([xA / 2. + xB / 2., yB / 2. + yA / 2.])

    n1 = np.array([yA - yB, xB - xA])  #normal vector to segA
    n2 = np.array([yC - yD, xD - xC])  #normal vector to segB

    test1 = np.dot(n2, O2A)
    test2 = np.dot(n2, O2B)

    if test1 * test2 > 0:
        bval = 0
        return bval

    test3 = np.dot(n1, O1C)
    test4 = np.dot(n1, O1D)

    if test3 * test4 > 0:
        bval = 0
        return bval

    #if colinear
    if test1 * test2 == 0 and test3 * test4 == 0 and np.linalg.det(np.hstack((n1.reshape((-1, )), n2.reshape(-1, )))) == 0:

        #projection on the axis O1O2
        O2O1 = np.array([xA / 2. + xB / 2., yB / 2. + yA / 2.]) - np.array([xD / 2. + xC / 2., yD / 2. + yC / 2.])
        O1A = np.dot(O2O1, (O2A - O2O1))
        O1B = np.dot(O2O1, (O2B - O2O1))
        O1C = np.dot(O2O1, O1C)
        O1D = np.dot(O2O1, O1D)

    #test if one point is included in the other segment (-> bval = 1)
        if (O1C - O1A) * (O1D - O1A) < 0:
            bval = 1
            return bval
        if (O1C - O1B) * (O1D - O1B) < 0:
            bval = 1
            return bval
        if (O1A - O1C) * (O1B - O1C) < 0:
            bval = 1
            return bval
        if (O1A - O1D) * (O1B - O1D) < 0:
            bval = 1
            return bval

    #test if the 2 segments have the same middle (-> bval = 1)
        if O2O1 == 0:
            bval = 1
            return bval

    #else
        bval = 0
        return bval

    return bval
