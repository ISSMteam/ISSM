import numpy as np

def outflow(md):
    """
    OUTFLOW - flag nodes on outflux boundary.
    
    Python translation of:
        function flag = outflow(md)
        A = md.mesh.segments(:,1);
        B = md.mesh.segments(:,2);
        Nx = -(md.mesh.y(A) - md.mesh.y(B));
        Ny =  md.mesh.x(A) - md.mesh.x(B);
        Vx = (md.initialization.vx(A) + md.initialization.vx(B))/2;
        Vy = (md.initialization.vy(A) + md.initialization.vy(B))/2;
        VdotN = Vx.*Nx + Vy.*Ny;
        flag(A(find(VdotN>0))) = 1;
    """

    # Extract segment endpoints (convert to 0-based indexing if needed)
    A = md.mesh.segments[:, 0].astype(int)
    B = md.mesh.segments[:, 1].astype(int)

    # Handle potential 1-based indexing from ISSM .mat imports
    if A.min() >= 1 and B.min() >= 1:
        A = A - 1
        B = B - 1

    # Compute outward normals (Nx, Ny)
    Nx = -(md.mesh.y[A] - md.mesh.y[B])
    Ny =  (md.mesh.x[A] - md.mesh.x[B])

    # Average velocities at segment endpoints
    Vx = 0.5 * (md.initialization.vx[A] + md.initialization.vx[B])
    Vy = 0.5 * (md.initialization.vy[A] + md.initialization.vy[B])

    # Dot product of velocity and outward normal
    VdotN = Vx * Nx + Vy * Ny

    # Initialize vertex flags
    flag = np.zeros(md.mesh.numberofvertices, dtype=bool)

    # Mark vertices A for which flow is outward (V·N > 0)
    out_idx = np.where(VdotN > 0)[0]
    flag[A[out_idx]] = True

    return flag
