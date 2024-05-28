function SetIceSheetBC(md) { 
//SETICESHEETBC - Create the boundary conditions for stressbalance and thermal models for an IceSheet with no Ice Front
//
//   Usage
//      md=SetIceSheetBC(md)
//
//   See also: SETICESHELFBC, SETMARINEICESHEETBC

	md.stressbalance.spcvx=NewArrayFill(md.mesh.numberofvertices,NaN); 
	md.stressbalance.spcvy=NewArrayFill(md.mesh.numberofvertices,NaN);
	md.stressbalance.spcvz=NewArrayFill(md.mesh.numberofvertices,NaN);
	md.stressbalance.referential=Create2DArray(md.mesh.numberofvertices,6);
	for(var i=0;i<md.mesh.numberofvertices;i++)for(var j=0;j<6;j++)md.stressbalance.referential[i][j]=NaN;
	md.stressbalance.loadingforce=Create2DArray(md.mesh.numberofvertices,3);
	for(var i=0;i<md.mesh.numberofvertices;i++)for(var j=0;j<3;j++)md.stressbalance.loadingforce[i][j]=0;

	//Node on dirichlet: 
	pos=ArrayFind(md.mesh.vertexonboundary,1);
	for(var i=0;i<pos.length;i++){
		md.stressbalance.spcvx[pos[i]]=0;
		md.stressbalance.spcvy[pos[i]]=0;
		md.stressbalance.spcvz[pos[i]]=0;
	}

	//Dirichlet values: 
	var obs=false;
	if((md.inversion.vx_obs.length == md.mesh.numberofvertices) & (md.inversion.vy_obs.length==md.mesh.numberofvertices))obs=true;

	if(obs==true){
		console.log('      boundary conditions for stressbalance model: setting spc as observed velocities');
		for(var i=0;i<pos.length;i++){
			md.stressbalance.spcvx[pos[i]]=md.inversion.vx_obs[pos[i]];
			md.stressbalance.spcvy[pos[i]]=md.inversion.vy_obs[pos[i]];
		}
	}
	else{
		console.log('      boundary conditions for stressbalance model: setting spc as zero');
	}

	//No ice front, do nothing. 
	
	//Initialize surface and basal forcings
	md.smb.initialize(md);
	md.basalforcings.initialize(md);

	//Deal with other boundary conditions
	if (isNaN(md.balancethickness.thickening_rate)){
		md.balancethickness.thickening_rate=NewArrayFill(md.mesh.numberofvertices,0);
		console.log('      no balancethickness.thickening_rate specified: values set as zero');
	}
		
	md.masstransport.spcthickness=NewArrayFill(md.mesh.numberofvertices,NaN);
	md.balancethickness.spcthickness=NewArrayFill(md.mesh.numberofvertices,NaN);
	md.damage.spcdamage=NewArrayFill(md.mesh.numberofvertices,NaN);

	if (md.initialization.temperature.length==md.mesh.numberofvertices){
		md.thermal.spctemperature=NewArrayFill(md.mesh.numberofvertices,NaN);
		if ('vertexonsurface' in md.mesh){
			pos=ArrayFind(md.mesh.vertexonsurface,1);
			for(var i=0;i<pos.length;i++)md.thermal.spctemperature[i]=md.initialization.temperature[i]; //impose observed temperature on surface
		}
		if (md.basalforcings.geothermalflux.length != md.mesh.numberofvertices){
			md.basalforcings.geothermalflux=NewArrayFill(md.mesh.numberofvertices,0);
		}
	}
	else{
		console.log('      no thermal boundary conditions created: no observed temperature found');
	}
}
