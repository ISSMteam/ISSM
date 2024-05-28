from MatlabFuncs import *
import numpy as np
from MeshProfileIntersection import *
from helpers import empty_nd_list


def process_qmu_response_data(md):
    '''
PROCESS_QMU_RESPONSE_DATA - process any data necessary for the solutions to process the data.

    Usage: md = process_qmu_response_data(md)

    See also PREQMU, PRESOLVE
'''

    # preliminary data
    process_mass_flux_profiles = 0
    num_mass_flux = 0

    # loop through response descriptors, and act accordingly
    for i in range(np.size(md.qmu.responsedescriptors)):

        # Do we have to process  mass flux profiles?
        if strncmpi(md.qmu.responsedescriptors[i], 'indexed_MassFlux', 16):
            num_mass_flux += 1
            process_mass_flux_profiles = 1

    # deal with mass flux profiles
    if process_mass_flux_profiles:
        # we need a profile of points on which to compute the mass_flux, is it here?
        if type(md.qmu.mass_flux_profiles) == float and np.isnan(md.qmu.mass_flux_profiles):
            raise RuntimeError('process_qmu_response_data error message: could not find a mass_flux exp profile!')

        if type(md.qmu.mass_flux_profiles) != list:
            raise RuntimeError('process_qmu_response_data error message: qmu_mass_flux_profiles field should be a list of domain outline names')

        if np.size(md.qmu.mass_flux_profiles) == 0:
            raise RuntimeError('process_qmu_response_data error message: qmu_mass_flux_profiles cannot be empty!')

        if num_mass_flux != np.size(md.qmu.mass_flux_profiles):
            raise RuntimeError('process_qmu_response_data error message: qmu_mass_flux_profiles should be of the same size as the number of MassFlux responses asked for in the Qmu analysis')

    # ok, process the domains named in qmu_mass_flux_profiles,
    #     to build a list of segments (MatArray)
        md.qmu.mass_flux_segments = empty_nd_list((num_mass_flux, 1))

        for i in range(num_mass_flux):
            md.qmu.mass_flux_segments[i] = np.array(MeshProfileIntersection(md.mesh.elements, md.mesh.x, md.mesh.y, md.qmu.mass_flux_profile_directory + '/' + md.qmu.mass_flux_profiles[i])[0])

    return md
