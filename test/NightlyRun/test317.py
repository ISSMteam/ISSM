#Test Name: SquareSheetConstrainedTranHO
from model import *
from socket import gethostname
from triangle import *
from setmask import *
from parameterize import *
from setflowequation import *
from solve import *


md = triangle(model(), '../Exp/Square.exp', 200000.)
md = setmask(md, '', '')
md = parameterize(md, '../Par/SquareSheetConstrained.py')
md.basalforcings.groundedice_melting_rate[:] = 5.
md.extrude(3, 1.)
md = setflowequation(md, 'HO', 'all')
md.cluster = generic('name', gethostname(), 'np', 3)
md.transient.requested_outputs = ['default', 'GroundedArea', 'FloatingArea', 'TotalFloatingBmb', 'TotalGroundedBmb']
md = solve(md, 'Transient')

#Fields and tolerances to track changes
field_names = ['Vx1', 'Vy1', 'Vz1', 'Vel1', 'Pressure1', 'Bed1', 'Surface1', 'Thickness1', 'Temperature1', 'BasalforcingsGroundediceMeltingRate1', 'GroundedArea1', 'FloatingArea1', 'TotalFloatingBmb1', 'TotalGroundedBmb1',
               'Vx2', 'Vy2', 'Vz2', 'Vel2', 'Pressure2', 'Bed2', 'Surface2', 'Thickness2', 'Temperature2', 'BasalforcingsGroundediceMeltingRate2', 'GroundedArea2', 'FloatingArea2', 'TotalFloatingBmb2', 'TotalGroundedBmb2',
               'Vx3', 'Vy3', 'Vz3', 'Vel3', 'Pressure3', 'Bed3', 'Surface3', 'Thickness3', 'Temperature3', 'BasalforcingsGroundediceMeltingRate3', 'GroundedArea3', 'FloatingArea3', 'TotalFloatingBmb2', 'TotalGroundedBmb2']
field_tolerances = [1e-09, 1e-09, 1e-09, 1e-09, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-12, 1e-12, 1e-12, 1e-12,
                    1e-09, 1e-09, 1e-09, 1e-09, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-12, 1e-12, 1e-12, 1e-12,
                    1e-09, 5e-10, 1e-09, 1e-09, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-12, 1e-12, 1e-12, 1e-12]
field_values = [md.results.TransientSolution[0].Vx,
                md.results.TransientSolution[0].Vy,
                md.results.TransientSolution[0].Vz,
                md.results.TransientSolution[0].Vel,
                md.results.TransientSolution[0].Pressure,
                md.results.TransientSolution[0].Base,
                md.results.TransientSolution[0].Surface,
                md.results.TransientSolution[0].Thickness,
                md.results.TransientSolution[0].Temperature,
                md.results.TransientSolution[0].BasalforcingsGroundediceMeltingRate,
                md.results.TransientSolution[0].GroundedArea,
                md.results.TransientSolution[0].FloatingArea,
                md.results.TransientSolution[0].TotalFloatingBmb,
                md.results.TransientSolution[0].TotalGroundedBmb,
                md.results.TransientSolution[1].Vx,
                md.results.TransientSolution[1].Vy,
                md.results.TransientSolution[1].Vz,
                md.results.TransientSolution[1].Vel,
                md.results.TransientSolution[1].Pressure,
                md.results.TransientSolution[1].Base,
                md.results.TransientSolution[1].Surface,
                md.results.TransientSolution[1].Thickness,
                md.results.TransientSolution[1].Temperature,
                md.results.TransientSolution[1].BasalforcingsGroundediceMeltingRate,
                md.results.TransientSolution[1].GroundedArea,
                md.results.TransientSolution[1].FloatingArea,
                md.results.TransientSolution[1].TotalFloatingBmb,
                md.results.TransientSolution[1].TotalGroundedBmb,
                md.results.TransientSolution[2].Vx,
                md.results.TransientSolution[2].Vy,
                md.results.TransientSolution[2].Vz,
                md.results.TransientSolution[2].Vel,
                md.results.TransientSolution[2].Pressure,
                md.results.TransientSolution[2].Base,
                md.results.TransientSolution[2].Surface,
                md.results.TransientSolution[2].Thickness,
                md.results.TransientSolution[2].Temperature,
                md.results.TransientSolution[2].BasalforcingsGroundediceMeltingRate,
                md.results.TransientSolution[2].GroundedArea,
                md.results.TransientSolution[2].FloatingArea,
                md.results.TransientSolution[2].TotalFloatingBmb,
                md.results.TransientSolution[2].TotalGroundedBmb]
