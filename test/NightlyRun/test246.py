#Test Name: SquareShelfTranSemic
import numpy as np
from model import *
from socket import gethostname
from triangle import triangle
from setmask import setmask
from parameterize import parameterize
from setflowequation import setflowequation
from solve import solve
from generic import generic
from SMBsemic import SMBsemic

md=triangle(model(),'../Exp/Square.exp',150000.)
md=setmask(md,'all','')
md=parameterize(md,'../Par/SquareShelf.py')

# Use of SMBsemic
md.smb = SMBsemic()
# initalize pdd fields
md.smb=md.smb.initialize(md)
md.smb.s0gcm=md.geometry.surface

# Okay, initialize attribute matrix
for field in ['dailytemperature','dailysnowfall','dailyrainfall','dailydsradiation','dailydlradiation','dailywindspeed','dailyairdensity','dailyairhumidity','dailypressure']:
    setattr(md.smb, field, np.zeros((md.mesh.numberofvertices+1,365+1)))

ONES=np.ones((md.mesh.numberofvertices,),dtype=float)
for iday in range(365+1):
	md.smb.dailytemperature[0:md.mesh.numberofvertices,iday]=252.8739*ONES
	md.smb.dailytemperature[md.mesh.numberofvertices,iday]=((iday)/12)
	md.smb.dailysnowfall[0:md.mesh.numberofvertices,iday]=8.5503e-09*ONES
	md.smb.dailysnowfall[md.mesh.numberofvertices,iday]=((iday)/12)
	md.smb.dailyrainfall[0:md.mesh.numberofvertices,iday]=1.7296e-09*ONES
	md.smb.dailyrainfall[md.mesh.numberofvertices,iday]=((iday)/12)
	md.smb.dailydsradiation[0:md.mesh.numberofvertices,iday]=128.1702*ONES
	md.smb.dailydsradiation[md.mesh.numberofvertices,iday]=((iday)/12)
	md.smb.dailydlradiation[0:md.mesh.numberofvertices,iday]=176.5667*ONES
	md.smb.dailydlradiation[md.mesh.numberofvertices,iday]=((iday)/12)
	md.smb.dailywindspeed[0:md.mesh.numberofvertices,iday]=6.0741*ONES
	md.smb.dailywindspeed[md.mesh.numberofvertices,iday]=((iday)/12)
	md.smb.dailyairdensity[0:md.mesh.numberofvertices,iday]=1.0729*ONES
	md.smb.dailyairdensity[md.mesh.numberofvertices,iday]=((iday)/12)
	md.smb.dailyairhumidity[0:md.mesh.numberofvertices,iday]=9.6667*10**(-4)*ONES
	md.smb.dailyairhumidity[md.mesh.numberofvertices,iday]=((iday)/12)
	md.smb.dailypressure[0:md.mesh.numberofvertices,iday]=7.7841*10**(4)*ONES
	md.smb.dailypressure[md.mesh.numberofvertices,iday]=((iday)/12)

# time steps and resolution
md.timestepping.time_step=0.5;
md.settings.output_frequency=1;
md.timestepping.final_time=1;

md.transient.issmb=1;
md.transient.ismasstransport=0;
md.transient.isstressbalance=0;
md.transient.isthermal=0;

md.transient.requested_outputs=['default','TemperatureSEMIC']
md.cluster=generic('name',oshostname(),'np',4) # 3 for the cluster
md=solve(md,'Transient')

#Fields and tolerances to track changes
field_names     =['TemperatureSEMIC1','SmbMassBalance1','TemperatureSEMIC2','SmbMassBalance2']
field_tolerances=[1e-13,1e-13,1e-13,1e-13]
field_values=[md.results.TransientSolution[0].TemperatureSEMIC,
			md.results.TransientSolution[0].SmbMassBalance,
			md.results.TransientSolution[1].TemperatureSEMIC,
			md.results.TransientSolution[1].SmbMassBalance,
			  ]
