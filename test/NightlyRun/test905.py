#Test Name: SquareNoDynHydrologyDCSmbCoupled
import numpy as np
from model import *
from socket import gethostname
from triangle import triangle
from setmask import setmask
from parameterize import parameterize
from solve import solve
from generic import generic
from SMBgradientscomponents import SMBgradientscomponents

md = triangle(model(), '../Exp/Square.exp', 100000.)
md = setmask(md, '', '')
md = parameterize(md, '../Par/SquareNoDyn.py')
md.cluster = generic('name', gethostname(), 'np', 1)

md.transient.ishydrology = True
md.transient.issmb = True
md.hydrology = hydrologydc()
md.hydrology = md.hydrology.initialize(md)
md.smb = SMBgradientscomponents()

md.hydrology.isefficientlayer = 1

md.hydrology.sedimentlimit_flag = 1
md.hydrology.sedimentlimit = 400.0
md.hydrology.sediment_transmitivity = 3.0 * np.ones((md.mesh.numberofvertices))
md.hydrology.mask_thawed_node = np.ones((md.mesh.numberofvertices))

md.hydrology.mask_eplactive_node = np.zeros((md.mesh.numberofvertices))
md.hydrology.epl_conductivity = 3.
md.hydrology.epl_initial_thickness = 20
md.hydrology.epl_colapse_thickness = 1.0e-3
md.hydrology.epl_thick_comp = 0
md.hydrology.epl_max_thickness = 1

md.hydrology.spcsediment_head = np.nan * np.ones((md.mesh.numberofvertices))
md.hydrology.spcepl_head = np.nan * np.ones((md.mesh.numberofvertices))

md.initialization.sediment_head = np.zeros((md.mesh.numberofvertices))
md.initialization.epl_head = np.zeros((md.mesh.numberofvertices))
md.initialization.epl_thickness = np.ones((md.mesh.numberofvertices))

#try to keep the different steps as common multipliers of eachother
#for the sake of precision the values of the model used as input should be on a shorter time-step-
#you can plot the results of this test and check how the time stepping is doing (see bellow)
md.hydrology.steps_per_step = 5
md.smb.steps_per_step = 10  #md.hydrology.steps_per_step
md.timestepping.time_step = 1.
md.timestepping.final_time = 20.0

smb_step = md.timestepping.time_step / md.smb.steps_per_step
hydro_step = md.timestepping.time_step / md.hydrology.steps_per_step
duration = np.arange(md.timestepping.start_time, md.timestepping.final_time + smb_step, smb_step)

ddf = 10.0e-3
md.smb.accuref = np.array([[0.5, 0.5], [md.timestepping.start_time, md.timestepping.final_time]])
md.smb.accualti = 0.0
md.smb.accugrad = np.array([[0., 0.], [md.timestepping.start_time, md.timestepping.final_time]])

#md.smb.runoffref = 9. * ddf * np.ones(np.shape(duration))  #constant input for testing
md.smb.runoffref = 0.9 * duration * ddf
md.smb.runoffref = np.vstack((md.smb.runoffref, duration))
md.smb.runoffalti = 0.0
md.smb.runoffgrad = np.array([[-6.5e-3 * ddf, -6.5e-3 * ddf], [md.timestepping.start_time, md.timestepping.final_time]])  #lapse rate *ddf*day per year

md.basalforcings.groundedice_melting_rate = np.zeros((md.mesh.numberofvertices))
md.basalforcings.floatingice_melting_rate = np.zeros((md.mesh.numberofvertices))

md = solve(md, 'Transient')
field_names = ['SedimentWaterHead1', 'EplWaterHead1', 'SedimentHeadResidual1',
               'SedimentWaterHead4', 'EplWaterHead4', 'SedimentHeadResidual4',
               'SedimentWaterHead5', 'EplWaterHead5', 'SedimentHeadResidual5',
               'SedimentWaterHead9', 'EplWaterHead9', 'SedimentHeadResidual9',
               'EplWaterHead20', 'EplWaterHeadSubstep20', 'SedimentWaterHead20',
               'SedimentWaterHeadSubstep20']
field_tolerances = [1e-13, 1e-13, 1e-13,
                    1e-13, 1e-13, 1e-13,
                    1e-13, 5e-12, 1e-11,
                    1e-13, 5e-12, 1e-11,
                    1e-13, 1e-13, 1e-13,
                    1e-13]
field_values = [md.results.TransientSolution[0].SedimentHead,
                md.results.TransientSolution[0].EplHead,
                md.results.TransientSolution[0].SedimentHeadResidual,
                md.results.TransientSolution[3].SedimentHead,
                md.results.TransientSolution[3].EplHead,
                md.results.TransientSolution[3].SedimentHeadResidual,
                md.results.TransientSolution[4].SedimentHead,
                md.results.TransientSolution[4].EplHead,
                md.results.TransientSolution[4].SedimentHeadResidual,
                md.results.TransientSolution[8].SedimentHead,
                md.results.TransientSolution[8].EplHead,
                md.results.TransientSolution[8].SedimentHeadResidual,
                md.results.TransientSolution[-1].EplHead,
                md.results.TransientSolution[-1].EplHeadSubstep,
                md.results.TransientSolution[-1].SedimentHead,
                md.results.TransientSolution[-1].SedimentHeadSubstep]


#===Stepping and plotting reference
# import matplotlib as mpl
# import matplotlib.pyplot as plt

# agregatedInput = np.zeros(len(md.results.TransientSolution))
# agregatedInputSub = np.zeros(len(md.results.TransientSolution))
# RefInput = np.zeros(len(md.results.TransientSolution))
# RefInputSub = np.zeros(len(md.results.TransientSolution))
# Input = np.zeros(len(md.results.TransientSolution))
# InputSub = np.zeros(len(md.results.TransientSolution))
# SedVol = np.zeros(len(md.results.TransientSolution))
# SedVolSub = np.zeros(len(md.results.TransientSolution))
# EplVol = np.zeros(len(md.results.TransientSolution))
# EplVolSub = np.zeros(len(md.results.TransientSolution))
# SmbTimer = np.zeros(len(md.results.TransientSolution))
# HydroTimer = np.zeros(len(md.results.TransientSolution))
# TimerSub = np.zeros(len(md.results.TransientSolution))
# HeadEVol = np.zeros(len(md.results.TransientSolution))
# sedstore = md.materials.rho_freshwater * md.constants.g * md.hydrology.sediment_porosity * md.hydrology.sediment_thickness * ((md.hydrology.sediment_compressibility / md.hydrology.sediment_porosity) + md.hydrology.water_compressibility)
# eplstore = md.materials.rho_freshwater * md.constants.g * md.hydrology.sediment_porosity * ((md.hydrology.sediment_compressibility / md.hydrology.sediment_porosity) + md.hydrology.water_compressibility)

# for i, res in enumerate(md.results.TransientSolution):

#     TimerSub[i] = res.time
#     RefInputSub[i] = max(0, res.time * 0.9 * ddf + (-6.5e-3 * 1000 * ddf))
#     try:
#         InputSub[i] = np.nanmean(res.SmbRunoffSubstep)
#         SedVolSub[i] = np.nanmean(res.SedimentHeadSubstep) * sedstore
#         EplVolSub[i] = np.nanmean(res.EplHeadSubstep) * eplstore * np.nanmean(res.HydrologydcEplThicknessSubstep)
#         substep = True
#         if i > 0:
#             agregatedInputSub[i:] += 0.5 * ((1 + smb_step) * np.nanmean(res.SmbRunoffSubstep) + (1 - smb_step) * np.nanmean(md.results.TransientSolution[i - 1].SmbRunoffSubstep)) * (TimerSub[i] - TimerSub[i - 1])
#         else:
#             agregatedInputSub[i:] += np.nanmean(res.SmbRunoffSubstep) * TimerSub[i]

#     except AttributeError:
#         substep = False
#     SmbTimer[i] = res.time - (0.5 * md.timestepping.time_step * (1 - smb_step))
#     HydroTimer[i] = res.time - (0.5 * md.timestepping.time_step * (1 - hydro_step))

#     RefInput[i] = max(0, (SmbTimer[i] * 0.9 * ddf + (-6.5e-3 * 1000 * ddf)))
#     Input[i] = np.nanmean(res.SmbRunoff)
#     SedVol[i] = np.nanmean(res.SedimentHead) * sedstore
#     EplVol[i] = np.nanmean(res.EplHead) * eplstore * np.nanmean(res.HydrologydcEplThickness)
#     if i > 0:
#         agregatedInput[i:] += 0.5 * ((1 + hydro_step) * np.nanmean(res.SmbRunoff) + (1 - hydro_step) * np.nanmean(md.results.TransientSolution[i - 1].SmbRunoff)) * (SmbTimer[i] - SmbTimer[i - 1])
#     else:
#         agregatedInput[i:] += np.nanmean(res.SmbRunoff) * SmbTimer[i]


# layout, ax = plt.subplots(2, 2, sharex=True, sharey=False, figsize=(10, 10))

# DiffTimer = (TimerSub[:-1] + 0.5 * (smb_step + hydro_step))
# Indiff = (DiffTimer) * 0.9 * ddf + (-6.5e-3 * 1000 * ddf)
# Indiff[np.where(Indiff < 0)] = 0

# ax[0, 0].plot(SmbTimer, RefInput - Input, 'g+')
# ax[0, 0].plot(DiffTimer, (Indiff - np.diff((SedVol+EplVol)) / md.timestepping.time_step), 'kx')

# ax[0, 1].plot(SmbTimer, RefInput, 'm+')
# ax[0, 1].plot(SmbTimer, Input, 'g+')
# ax[0, 1].plot(DiffTimer, np.diff(SedVol+EplVol) / md.timestepping.time_step, 'r')

# ax[1, 0].plot(HydroTimer, (agregatedInput - SedVol - EplVol), 'b+')

# ax[1, 1].plot(HydroTimer, SedVol + EplVol, 'r+')
# ax[1, 1].plot(SmbTimer, agregatedInput, 'b+')

# if substep:
#     SubDiffTimer = (SmbTimer[1:])
#     SubIndiff = (SubDiffTimer) * 0.9 * ddf + (-6.5e-3 * 1000 * ddf)
#     SubIndiff[np.where(SubIndiff < 0)] = 0

#     ax[0, 0].plot(TimerSub, RefInputSub - InputSub, 'm')
#     ax[0, 0].plot(SubDiffTimer, (SubIndiff - np.diff(SedVolSub + EplVolSub) / md.timestepping.time_step), 'y')

#     ax[0, 1].plot(TimerSub, RefInputSub, 'm')
#     ax[0, 1].plot(TimerSub, InputSub, 'g')

#     ax[1, 0].plot(TimerSub, (agregatedInputSub - SedVolSub - EplVolSub), 'r')

#     ax[1, 1].plot(TimerSub, SedVolSub + EplVolSub, 'r')
#     ax[1, 1].plot(TimerSub, SedVolSub, 'r', alpha=0.5)
#     ax[1, 1].plot(TimerSub, EplVolSub, 'r', alpha=0.5)
#     ax[1, 1].plot(TimerSub, agregatedInputSub, 'b')
