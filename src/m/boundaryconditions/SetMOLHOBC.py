import numpy as np

def SetMOLHOBC(md):
    """
    SETMOLHOBC - Create the boundary conditions for stressbalance for MOLHO: VxBase, VyBase, VxShear, VyShear

       Usage:
          md = SetIceShelfBC(md, varargin)

       Example:
          md = SetIceShelfBC(md)

    """

    #node on Dirichlet (boundary and ~icefront)
    md.stressbalance.spcvx_base = md.stressbalance.spcvx
    md.stressbalance.spcvy_base = md.stressbalance.spcvy
    md.stressbalance.spcvx_shear = np.nan * md.stressbalance.spcvx
    md.stressbalance.spcvy_shear = np.nan * md.stressbalance.spcvy

    return md
