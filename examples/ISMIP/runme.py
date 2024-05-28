import numpy as np
from model import *
from squaremesh import squaremesh
from plotmodel import plotmodel
from export_netCDF import export_netCDF
from loadmodel import loadmodel
from setmask import setmask
from parameterize import parameterize
from setflowequation import setflowequation
from socket import gethostname
from solve import solve
#which steps to perform; steps are from 1 to 8
#step 7 is specific to ISMIPA
#step 8 is specific to ISMIPF

steps = [1]

# parameter file to be used, choose between IsmipA.py or IsmipF.py
ParamFile = 'IsmipA.py'

#Run Steps

#Mesh Generation #1
if 1 in steps:
    print("Now generating the mesh")
    #initialize md as a new model #help(model)
    #->

    # generate a squaremesh #help(squaremesh)
    # Side is 80 km long with 20 points
    #->

    # plot the given mesh #plotdoc
    #->

    # save the given model
    #->

#Masks #2
if 2 in steps:
    print("Setting the masks")
    # load the preceding step #help(loadmodel)
    # path is given by the organizer with the name of the given step
    #->

    # set the mask #help(setmask)
    # all MISMIP nodes are grounded
    #->

    # plot the given mask #md.mask to locate the field
    #->

    # save the given model
    #->

#Parameterization #3
if 3 in steps:
    print("Parameterizing")
    # load the preceding step #help(loadmodel)
    # path is given by the organizer with the name of the given step
    #->

    # parametrize the model # help(parameterize)
    # you will need to fil-up the parameter file (given by the
    # ParamFile variable)
    #->

    # save the given model
    #->


#Extrusion #4
if 4 in steps:
    print("Extruding")
    # load the preceding step #help(loadmodel)
    # path is given by the organizer with the name of the given step
    #->

    # vertically extrude the preceding mesh #help(extrude)
    # only 5 layers exponent 1
    #->

    # plot the 3D geometry #plotdoc
    #->

    # save the given model
    #->

#Set the flow computing method #5
if 5 in steps:
    print("setting flow approximation")
    # load the preceding step #help(loadmodel)
    # path is given by the organizer with the name of the given step
    #->

    # set the approximation for the flow computation #help(setflowequation)
    # We will be using the Higher Order Model (HO)
    #->

    # save the given model
    #->

#Set Boundary Conditions #6
if 6 in steps:
    print("setting boundary conditions")
    # load the preceding step #help(loadmodel)
    # path is given by the organizer with the name of the given step
    #->

    # dirichlet boundary condition are known as SPCs
    # ice frozen to the base, no velocity   #md.stressbalance
    # SPCs are initialized at NaN one value per vertex
    #->

    #->

    #->

    # extract the nodenumbers at the base #md.mesh.vertexonbase
    #->

    # set the sliding to zero on the bed (Vx and Vy)
    #->

    #->

    # periodic boundaries have to be fixed on the sides
    # Find the indices of the sides of the domain, for x and then for y
    # for x
    # create maxX, list of indices where x is equal to max of x
    #->

    # create minX, list of indices where x is equal to min of x
    #->

    # for y
    # create maxY, list of indices where y is equal to max of y
    # but not where x is equal to max or min of x
    # (i.e, indices in maxX and minX should be excluded from maxY and minY)
    #->

    # create minY, list of indices where y is equal to max of y
    # but not where x is equal to max or min of x
    #->

    # set the node that should be paired together, minX with maxX and minY with maxY
    # #md.stressbalance.vertex_pairing
    #->

    if ParamFile=='IsmipF.py':
        # if we are dealing with IsmipF the solution is in masstransport
        md.masstransport.vertex_pairing = md.stressbalance.vertex_pairing

    # save the given model
    #->

#Solving #7
if 7 in steps:
    print("running the solver for the A case")
    # load the preceding step #help(loadmodel)
    # path is given by the organizer with the name of the given step
    #->

    # Set cluster #md.cluster
    # generic parameters #help(generic)
    # set only the name and number of process
    #->

    # Set which control message you want to see #help(verbose)
    #->

    # Solve #help(solve)
    # we are solving a StressBalanc
    #->

    # save the given model
    #->

    # plot the surface velocities #plotdoc
    #->

#Solving #8
if 8 in steps:
    print("running the solver for the F case")
    # load the preceding step #help(loadmodel)
    # path is given by the organizer with the name of the given step
    #->

    # Set cluster #md.cluster
    # generic parameters #help(generic)
    # set only the name and number of process
    #->

    # Set which control message you want to see #help(verbose)
    #->

    # set the transient model to ignore the thermal model
    # #md.transient
    #->

    # define the timestepping scheme
    # everything here should be provided in years #md.timestepping
    # give the length of the time_step (4 years)
    #->

    # give final_time (20*4 years time_steps)
    #->

    # Solve #help(solve)
    # we are solving a TransientSolution
    #->

    # save the given model
    #->

    # plot the surface velocities #plotdoc
    #->
