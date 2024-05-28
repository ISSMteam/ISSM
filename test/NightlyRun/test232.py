#Test Name: SquareShelfTherTranForcTemp
import numpy as np
from model import *
from socket import gethostname
from triangle import *
from setmask import *
from parameterize import *
from setflowequation import *
from solve import *


md = triangle(model(), '../Exp/Square.exp', 180000.)
md = setmask(md, 'all', '')
md = parameterize(md, '../Par/SquareShelf.py')
md.extrude(3, 1.)
md = setflowequation(md, 'SSA', 'all')
md.cluster = generic('name', gethostname(), 'np', 3)
md.thermal.spctemperature = np.vstack((np.vstack((md.thermal.spctemperature,
                                                  md.thermal.spctemperature + 5.,
                                                  md.thermal.spctemperature + 10.,
                                                  md.thermal.spctemperature + 15.)).T, [1.5, 2.5, 3.5, 4.]))
md.timestepping.time_step = 1.
md.timestepping.final_time = 4.
md.transient.isstressbalance = False
md.transient.ismasstransport = False
md.transient.issmb = True
md.transient.isthermal = True
md.transient.isgroundingline = False
md = solve(md, 'Transient')

#Fields and tolerances to track changes
field_names = ['Temperature1', 'BasalforcingsGroundediceMeltingRate1',
               'Temperature2', 'BasalforcingsGroundediceMeltingRate2',
               'Temperature3', 'BasalforcingsGroundediceMeltingRate3',
               'Temperature4', 'BasalforcingsGroundediceMeltingRate4']
field_tolerances = [1e-13, 1e-6, 1e-13, 1e-6, 1e-13, 1e-6, 1e-13, 1e-6]
field_values = [md.results.TransientSolution[0].Temperature,
                md.results.TransientSolution[0].BasalforcingsGroundediceMeltingRate,
                md.results.TransientSolution[1].Temperature,
                md.results.TransientSolution[1].BasalforcingsGroundediceMeltingRate,
                md.results.TransientSolution[2].Temperature,
                md.results.TransientSolution[2].BasalforcingsGroundediceMeltingRate,
                md.results.TransientSolution[3].Temperature,
                md.results.TransientSolution[3].BasalforcingsGroundediceMeltingRate]
