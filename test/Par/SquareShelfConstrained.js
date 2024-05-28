function parameterize(md){

	//Geometry
	var hmin=300;
	var hmax=1000;
	var ymin=ArrayMin(md.mesh.y);
	var ymax=ArrayMax(md.mesh.y);
	var xmin=ArrayMin(md.mesh.x);
	var xmax=ArrayMax(md.mesh.x);
	
	md.geometry.thickness=NewArrayFill(md.mesh.numberofvertices,0);
	md.geometry.base=NewArrayFill(md.mesh.numberofvertices,0);
	md.geometry.surface=NewArrayFill(md.mesh.numberofvertices,0);
	md.geometry.bed=NewArrayFill(md.mesh.numberofvertices,0);

	for(i=0;i<md.mesh.numberofvertices;i++){
		md.geometry.thickness[i]=hmax+(hmin-hmax)*(md.mesh.y[i]-ymin)/(ymax-ymin)+0.1*(hmin-hmax)*(md.mesh.x[i]-xmin)/(xmax-xmin);
		md.geometry.base[i]=-md.materials.rho_ice/md.materials.rho_water*md.geometry.thickness[i];
		md.geometry.surface[i]=md.geometry.base[i]+md.geometry.thickness[i];
		md.geometry.bed[i]=md.geometry.base[i]-10;
	}
	
	//Initial velocity: no ncreader for now, so we just load arrays.
	/*x     = transpose(ncread('../Data/SquareShelfConstrained.nc','x'));
	y     = transpose(ncread('../Data/SquareShelfConstrained.nc','y'));
	vx    = transpose(ncread('../Data/SquareShelfConstrained.nc','vx'));
	vy    = transpose(ncread('../Data/SquareShelfConstrained.nc','vy'));
	index = transpose(ncread('../Data/SquareShelfConstrained.nc','index'));*/
	
	md.initialization.vx=InterpFromMeshToMesh2d(index,x,y,vx,md.mesh.x,md.mesh.y);
	md.initialization.vy=InterpFromMeshToMesh2d(index,x,y,vy,md.mesh.x,md.mesh.y);
	md.initialization.vel=ArrayMag(md.initialization.vx,md.initialization.vy);
	md.initialization.vz=NewArrayFill(md.mesh.numberofvertices,0);
	md.initialization.pressure=NewArrayFill(md.mesh.numberofvertices,0);

	//Materials
	md.initialization.temperature=NewArrayFill(md.mesh.numberofvertices,273-20);
	md.materials.rheology_B=paterson(md.initialization.temperature);
	md.materials.rheology_n=NewArrayFill(md.mesh.numberofelements,3);

	//Surface mass balance and basal melting
	md.smb.mass_balance=NewArrayFill(md.mesh.numberofvertices,10);
	md.basalforcings.floatingice_melting_rate=NewArrayFill(md.mesh.numberofvertices,5.0);
	md.basalforcings.groundedice_melting_rate=NewArrayFill(md.mesh.numberofvertices,5.0);

	//Friction
	md.friction.coefficient=NewArrayFill(md.mesh.numberofvertices,20);
	for(var i=0;i<md.mesh.numberofvertices;i++)if(md.mask.groundedice_levelset[i]<0)md.friction.coefficient[i]=0;
	md.friction.p=NewArrayFill(md.mesh.numberofelements,1);
	md.friction.q=NewArrayFill(md.mesh.numberofelements,1);

	//Numerical parameters
	md.masstransport.stabilization=1;
	md.thermal.stabilization=1;
	md.verbose=new verbose(0);
	md.settings.waitonlock=30;
	md.stressbalance.restol=0.05;
	md.stressbalance.reltol=0.05;
	md.steadystate.reltol=0.05;
	md.stressbalance.abstol=NaN;
	md.timestepping.time_step=1;
	md.timestepping.final_time=3;

	//Deal with boundary conditions:
	SetIceShelfBC(md);

	//Change name so that no tests have the same name
	md.miscellaneous.name='test101';
}
