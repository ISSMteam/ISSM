def cfl_step(md, vx, vy):
    """cfl_step - return the maximum time step for the model in years

    Dt < 0.5 / (u/Dx + v/Dy)

    Usage:
        maxtime = cfl_step(md, vx, vy)

    Example:
        dt=cfl_step(md,md.results.StressbalanceSolution.Vx,md.results.StressbalanceSolution.Vy)
    """

    import numpy as np
    # Return the maximum time step for the model in years
    if vx.shape[0] != md.mesh.numberofvertices and vy.shape[0] != md.mesh.numberofvertices:
        sys.exit('timesteps error message: size of velocity components must be the same as md.mesh.numberofvertices')


    index = md.mesh.elements
    edgex = np.nanmax(md.mesh.x[index-1], axis=1) - np.nanmin(md.mesh.x[index-1], axis=1)
    edgey = np.nanmax(md.mesh.y[index-1], axis=1) - np.nanmin(md.mesh.y[index-1], axis=1)
    vx = np.nanmax(np.abs(vx[index-1]), axis=1)
    vy = np.nanmax(np.abs(vy[index-1]), axis=1)
    maxtime = 1/2 * np.nanmin(1/(vx/edgex + vy/edgey))
    return maxtime
