import numpy as np
from paterson import paterson
from netCDF4 import Dataset
from InterpFromGridToMesh import InterpFromGridToMesh

#Name and Coordinate system
md.miscellaneous.name = 'SeaRISEgreenland'
md.mesh.epsg = 3413

print('   Loading SeaRISE data from NetCDF')
ncdata = Dataset('../Data/Greenland_5km_dev1.2.nc', mode='r')
x1 = np.squeeze(ncdata.variables['x1'][:].data)
y1 = np.squeeze(ncdata.variables['y1'][:].data)
usrf = np.squeeze(ncdata.variables['usrf'][:].data)
topg = np.squeeze(ncdata.variables['topg'][:].data)
velx = np.squeeze(ncdata.variables['surfvelx'][:].data)
vely = np.squeeze(ncdata.variables['surfvely'][:].data)
temp = np.squeeze(ncdata.variables['airtemp2m'][:].data)
smb = np.squeeze(ncdata.variables['smb'][:].data)
gflux = np.squeeze(ncdata.variables['bheatflx'][:].data)
ncdata.close()

print('   Interpolating surface and bedrock')
md.geometry.base = InterpFromGridToMesh(x1, y1, topg, md.mesh.x, md.mesh.y, 0)
md.geometry.surface = InterpFromGridToMesh(x1, y1, usrf, md.mesh.x, md.mesh.y, 0)

print('   Constructing thickness')
md.geometry.thickness = md.geometry.surface - md.geometry.base

#Set min thickness to 1 meter
pos0 = np.nonzero(md.geometry.thickness <= 0)
md.geometry.thickness[pos0] = 1
md.geometry.surface = md.geometry.thickness + md.geometry.base

print('   Interpolating velocities ')
md.inversion.vx_obs = InterpFromGridToMesh(x1, y1, velx, md.mesh.x, md.mesh.y, 0)
md.inversion.vy_obs = InterpFromGridToMesh(x1, y1, vely, md.mesh.x, md.mesh.y, 0)
md.inversion.vel_obs = np.sqrt(md.inversion.vx_obs**2 + md.inversion.vy_obs**2)
md.initialization.vx = md.inversion.vx_obs
md.initialization.vy = md.inversion.vy_obs
md.initialization.vz = np.zeros((md.mesh.numberofvertices))
md.initialization.vel = md.inversion.vel_obs

print('   Interpolating temperatures')
md.initialization.temperature = InterpFromGridToMesh(x1, y1, temp, md.mesh.x, md.mesh.y, 0) + 273.15

print('   Interpolating surface mass balance')
md.smb.mass_balance = InterpFromGridToMesh(x1, y1, smb, md.mesh.x, md.mesh.y, 0)
md.smb.mass_balance = md.smb.mass_balance * md.materials.rho_water / md.materials.rho_ice

print('   Construct basal friction parameters')
md.friction.coefficient = 30 * np.ones((md.mesh.numberofvertices))
pos = np.nonzero(md.mask.ocean_levelset < 0)
md.friction.coefficient[pos] = 0  #no friction applied on floating ice
md.friction.p = np.ones((md.mesh.numberofelements))
md.friction.q = np.ones((md.mesh.numberofelements))

print('   Construct ice rheological properties')
md.materials.rheology_n = 3 * np.ones((md.mesh.numberofelements))
md.materials.rheology_B = paterson(md.initialization.temperature)
md.friction.q = np.ones((md.mesh.numberofelements))
md.friction.p = np.ones((md.mesh.numberofelements))

print('   Set other boundary conditions')
md.mask.ice_levelset[np.nonzero(md.mesh.vertexonboundary == 1)] = 0
md.basalforcings.floatingice_melting_rate = np.zeros((md.mesh.numberofvertices))
md.basalforcings.groundedice_melting_rate = np.zeros((md.mesh.numberofvertices))
#impose observed temperature on surface
md.thermal.spctemperature = md.initialization.temperature
md.masstransport.spcthickness = np.nan * np.ones((md.mesh.numberofvertices))

print('   Set geothermal heat flux')
md.basalforcings.geothermalflux = InterpFromGridToMesh(x1, y1, gflux, md.mesh.x, md.mesh.y, 0)

print('   Set Pressure')
md.initialization.pressure = md.materials.rho_ice * md.constants.g * md.geometry.thickness

print('   Single point constraints')
#Initialize single point constraint arrays
md.stressbalance.referential = np.nan * np.ones((md.mesh.numberofvertices, 6))
md.stressbalance.spcvx = np.nan * np.ones((md.mesh.numberofvertices))
md.stressbalance.spcvy = np.nan * np.ones((md.mesh.numberofvertices))
md.stressbalance.spcvz = np.nan * np.ones((md.mesh.numberofvertices))
