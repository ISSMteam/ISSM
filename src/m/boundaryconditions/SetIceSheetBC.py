import numpy as np


def SetIceSheetBC(md):
    """SETICESHEETBC - Create the boundary conditions for stressbalance and thermal models for an IceSheet with no Ice Front

    Usage:
        md = SetIceSheetBC(md)

    See also: SETICESHELFBC, SETMARINEICESHEETBC
    """

    # node on Dirichlet
    pos = np.nonzero(md.mesh.vertexonboundary)
    md.stressbalance.spcvx = np.nan * np.ones((md.mesh.numberofvertices))
    md.stressbalance.spcvy = np.nan * np.ones((md.mesh.numberofvertices))
    md.stressbalance.spcvz = np.nan * np.ones((md.mesh.numberofvertices))
    md.stressbalance.spcvx[pos] = 0
    md.stressbalance.spcvy[pos] = 0
    md.stressbalance.spcvz[pos] = 0
    md.stressbalance.referential = np.nan * np.ones((md.mesh.numberofvertices, 6))
    md.stressbalance.loadingforce = 0 * np.ones((md.mesh.numberofvertices, 3))

    # Dirichlet Values
    if isinstance(md.inversion.vx_obs, np.ndarray) and np.size(md.inversion.vx_obs, axis=0) == md.mesh.numberofvertices and isinstance(md.inversion.vy_obs, np.ndarray) and np.size(md.inversion.vy_obs, axis=0) == md.mesh.numberofvertices:
        print('      boundary conditions for stressbalance model: spc set as observed velocities')
        md.stressbalance.spcvx[pos] = md.inversion.vx_obs[pos]
        md.stressbalance.spcvy[pos] = md.inversion.vy_obs[pos]
    else:
        print('      boundary conditions for stressbalance model: spc set as zero')

    # No ice front -> do nothing

    # Initialize surface and basal forcings
    md.smb.initialize(md)
    md.basalforcings.initialize(md)

    # Initialize ocean forcings and sealevel
    md.dsl.initialize(md)

    # Deal with other boundary conditions
    if np.all(np.isnan(md.balancethickness.thickening_rate)):
        md.balancethickness.thickening_rate = np.zeros((md.mesh.numberofvertices))
        print('      no balancethickness.thickening_rate specified: values set as zero')
    md.masstransport.spcthickness = np.nan * np.ones((md.mesh.numberofvertices))
    md.balancethickness.spcthickness = np.nan * np.ones((md.mesh.numberofvertices))
    md.damage.spcdamage = np.nan * np.ones((md.mesh.numberofvertices))

    if isinstance(md.initialization.temperature, np.ndarray) and np.size(md.initialization.temperature, axis=0) == md.mesh.numberofvertices:
        md.thermal.spctemperature = np.nan * np.ones((md.mesh.numberofvertices))
        if hasattr(md.mesh, 'vertexonsurface'):
            pos = np.nonzero(md.mesh.vertexonsurface)[0]
            md.thermal.spctemperature[pos] = md.initialization.temperature[pos]  # impose observed temperature on surface
        if not isinstance(md.basalforcings.geothermalflux, np.ndarray) or not np.size(md.basalforcings.geothermalflux) == md.mesh.numberofvertices:
            md.basalforcings.geothermalflux = 50.0 * pow(10, -3) * np.ones((md.mesh.numberofvertices))  # 50 mW/m^2
    else:
        print('      no thermal boundary conditions created: no observed temperature found')

    return md
