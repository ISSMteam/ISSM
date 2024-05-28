#Test Name: SquareShelfConstrainedStressMOLHO2d
from model import *
from socket import gethostname
from triangle import triangle
from setmask import setmask
from parameterize import parameterize
from setflowequation import setflowequation
from solve import solve
from massfluxatgate import massfluxatgate
from generic import generic
from SetMOLHOBC import SetMOLHOBC

md = triangle(model(), '../Exp/Square.exp', 50000)
md = setmask(md, 'all', '')
md = parameterize(md, '../Par/SquareShelfConstrained.py')
md = setflowequation(md, 'MOLHO', 'all')
md.cluster = generic('name', gethostname(), 'np', 2)
#outputs
#FIXME compute the stress components for MOLHO
md.stressbalance.requested_outputs = ['default', 'VxSurface', 'VySurface', 'VxShear', 'VyShear', 'VxBase', 'VyBase', 'MassFlux1', 'MassFlux2', 'MassFlux3', 'MassFlux4', 'MassFlux5', 'MassFlux6']
#md.stressbalance.requested_outputs = ['default', 'DeviatoricStressxx', 'DeviatoricStressyy', 'DeviatoricStressxy', 'MassFlux1', 'MassFlux2', 'MassFlux3', 'MassFlux4', 'MassFlux5', 'MassFlux6']
md.outputdefinition.definitions = [massfluxatgate('name', 'MassFlux1', 'profilename', '../Exp/MassFlux1.exp', 'definitionstring', 'Outputdefinition1'),
                                   massfluxatgate('name', 'MassFlux2', 'profilename', '../Exp/MassFlux2.exp', 'definitionstring', 'Outputdefinition2'),
                                   massfluxatgate('name', 'MassFlux3', 'profilename', '../Exp/MassFlux3.exp', 'definitionstring', 'Outputdefinition3'),
                                   massfluxatgate('name', 'MassFlux4', 'profilename', '../Exp/MassFlux4.exp', 'definitionstring', 'Outputdefinition4'),
                                   massfluxatgate('name', 'MassFlux5', 'profilename', '../Exp/MassFlux5.exp', 'definitionstring', 'Outputdefinition5'),
                                   massfluxatgate('name', 'MassFlux6', 'profilename', '../Exp/MassFlux6.exp', 'definitionstring', 'Outputdefinition6')]
md = SetMOLHOBC(md)
md = solve(md, 'Stressbalance')

#Fields and tolerances to track changes
field_names = ['Vx', 'Vy', 'Vel', 'Pressure', 'VxSurface', 'VySurface', 'VxShear', 'VyShear', 'VxBase', 'VyBase',
               'MassFlux1', 'MassFlux2', 'MassFlux3', 'MassFlux4', 'MassFlux5', 'MassFlux6']
#'DeviatoricStressxx', 'DeviatoricStressyy', 'DeviatoricStressxy',
field_tolerances = [5e-13, 6e-13, 6e-13, 1e-13, 5e-13, 6e-13, 2e-13, 4e-13, 5e-13, 6e-13,
                    2e-13, 1e-13, 2e-13,
                    1e-13, 1e-13, 1e-13]
                    #1e-13, 1e-13, 1e-13]
field_values = [md.results.StressbalanceSolution.Vx,
                md.results.StressbalanceSolution.Vy,
                md.results.StressbalanceSolution.Vel,
                md.results.StressbalanceSolution.Pressure,
                md.results.StressbalanceSolution.VxSurface,
                md.results.StressbalanceSolution.VySurface,
                md.results.StressbalanceSolution.VxShear,
                md.results.StressbalanceSolution.VyShear,
                md.results.StressbalanceSolution.VxBase,
                md.results.StressbalanceSolution.VyBase,
                md.results.StressbalanceSolution.MassFlux1,
                md.results.StressbalanceSolution.MassFlux2,
                md.results.StressbalanceSolution.MassFlux3,
                md.results.StressbalanceSolution.MassFlux4,
                md.results.StressbalanceSolution.MassFlux5,
                md.results.StressbalanceSolution.MassFlux6]
#md.results.StressbalanceSolution.DeviatoricStressxx,
#md.results.StressbalanceSolution.DeviatoricStressyy,
#md.results.StressbalanceSolution.DeviatoricStressxy,
