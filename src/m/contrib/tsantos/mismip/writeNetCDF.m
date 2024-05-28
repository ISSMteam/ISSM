function writeNetCDF(md0,md,step,dt,ncfile),

 	xGL={};
   yGL={};
   iceThicknessGL={};
   uBaseGL={};
   vBaseGL={};
   iceVolume=[];
   iceVAF=[];
   groundedArea=[];
   time=[];

	%Inserting time 0. md0 must be last experiment (e.g., Ice1r for Ice1ra)
	if(isfield(md0.results.TransientSolution,'MeshElements'))
      x     = md0.results.TransientSolution(end).MeshX;
      y     = md0.results.TransientSolution(end).MeshY;
   else
      x     = md0.mesh.x;
      y     = md0.mesh.y;
   end	
	time(1)					= 0;
	[xgl_step ygl_step]	= gl_position(md0,length(md0.results.TransientSolution),0);
	xGL{1}					= xgl_step;
	yGL{1}					= ygl_step;
	iceVolume(1)			= md0.results.TransientSolution(end).IceVolume;
	iceVAF(1)				= md0.results.TransientSolution(end).IceVolumeAboveFloatation;
	groundedArea(1)		= md0.results.TransientSolution(end).GroundedArea;
	iceThicknessGL{1}		= griddata(x,y,md0.results.TransientSolution(end).Thickness,xgl_step,ygl_step);
	uBaseGL{1}				= griddata(x,y,md0.results.TransientSolution(end).Vx,xgl_step,ygl_step);
	vBaseGL{1}				= griddata(x,y,md0.results.TransientSolution(end).Vy,xgl_step,ygl_step);

   for i=2:length(step),
		if(isfield(md.results.TransientSolution,'MeshElements'))
			x     = md.results.TransientSolution(step(i)).MeshX;
			y     = md.results.TransientSolution(step(i)).MeshY;
		else
			x     = md.mesh.x;
			y     = md.mesh.y;
		end	
		time(i)=md.results.TransientSolution(step(i)).time;	
		[xgl_step ygl_step]=gl_position(md,step(i),0);
      xGL{i}=xgl_step;
      yGL{i}=ygl_step;
      iceVolume(i)=md.results.TransientSolution(step(i)).IceVolume;
      iceVAF(i)=md.results.TransientSolution(step(i)).IceVolumeAboveFloatation;
      groundedArea(i)=md.results.TransientSolution(step(i)).GroundedArea;
      iceThicknessGL{i}=griddata(x,y,md.results.TransientSolution(step(i)).Thickness,xgl_step,ygl_step);
      uBaseGL{i}=griddata(x,y,md.results.TransientSolution(step(i)).Vx,xgl_step,ygl_step);
      vBaseGL{i}=griddata(x,y,md.results.TransientSolution(step(i)).Vy,xgl_step,ygl_step);
   end
   uSurfaceGL=uBaseGL;
   vSurfaceGL=vBaseGL;
   uMeanGL=uBaseGL;
   vMeanGL=vBaseGL;
	
	 %Create netcdf
   mode = netcdf.getConstant('NETCDF4');
   mode = bitor(mode,netcdf.getConstant('CLASSIC_MODEL'));
   ncid=netcdf.create(ncfile,mode);

   %General attributes
   netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Author','Helene Seroussi (helene.seroussi@jpl.nasa.gov)');
   netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Model','ISSM (Ice Sheet System Model)');
   netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Notes','Experiments performed at Caltech Jet Propulsion Laboratory, Pasadena');
   netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Date',date());

   %Define dimensions
   gl_id    = netcdf.defDim(ncid,'nPointGL',netcdf.getConstant('NC_UNLIMITED'));
   ntime_id  = netcdf.defDim(ncid,'nTime',length(time));

   %Define variables
   volume_var_id = netcdf.defVar(ncid,'iceVolume','NC_FLOAT',[ntime_id]);
   VAF_var_id = netcdf.defVar(ncid,'iceVAF','NC_FLOAT',[ntime_id]);
   grounded_var_id = netcdf.defVar(ncid,'groundedArea','NC_FLOAT',[ntime_id]);
   xgl_var_id = netcdf.defVar(ncid,'xGL','NC_FLOAT',[ntime_id,gl_id]);
   ygl_var_id = netcdf.defVar(ncid,'yGL','NC_FLOAT',[ntime_id,gl_id]);
   thickness_var_id = netcdf.defVar(ncid,'iceThicknessGL','NC_FLOAT',[ntime_id,gl_id]);
   ubase_var_id = netcdf.defVar(ncid,'uBaseGL','NC_FLOAT',[ntime_id,gl_id]);
   vbase_var_id = netcdf.defVar(ncid,'vBaseGL','NC_FLOAT',[ntime_id,gl_id]);
   usurface_var_id = netcdf.defVar(ncid,'uSurfaceGL','NC_FLOAT',[ntime_id,gl_id]);
   vsurface_var_id = netcdf.defVar(ncid,'vSurfaceGL','NC_FLOAT',[ntime_id,gl_id]);
   umean_var_id = netcdf.defVar(ncid,'uMeanGL','NC_FLOAT',[ntime_id,gl_id]);
   vmean_var_id = netcdf.defVar(ncid,'vMeanGL','NC_FLOAT',[ntime_id,gl_id]);
   time_var_id = netcdf.defVar(ncid,'time','NC_FLOAT',[ntime_id]);

	%Define default fill values
   if(false),
		netcdf.defVarFill(ncid,xgl_var_id,false,NaN);
		netcdf.defVarFill(ncid,ygl_var_id,false,NaN);
		netcdf.defVarFill(ncid,thickness_var_id,false,NaN);
		netcdf.defVarFill(ncid,ubase_var_id,false,NaN);
		netcdf.defVarFill(ncid,vbase_var_id,false,NaN);
		netcdf.defVarFill(ncid,usurface_var_id,false,NaN);
		netcdf.defVarFill(ncid,vsurface_var_id,false,NaN);
		netcdf.defVarFill(ncid,umean_var_id,false,NaN);
		netcdf.defVarFill(ncid,vmean_var_id,false,NaN);
	end
   
	netcdf.endDef(ncid);
	
	%Write variables
   netcdf.putVar(ncid,volume_var_id,iceVolume);
   netcdf.putVar(ncid,VAF_var_id,iceVAF);
   netcdf.putVar(ncid,grounded_var_id,groundedArea);
   for i=1:length(time),
      netcdf.putVar(ncid,xgl_var_id,[i-1,0],[1,length(xGL{i})],xGL{i}');
      netcdf.putVar(ncid,ygl_var_id,[i-1,0],[1,length(xGL{i})],yGL{i}');
      netcdf.putVar(ncid,thickness_var_id,[i-1,0],[1,length(xGL{i})],iceThicknessGL{i}');
      netcdf.putVar(ncid,ubase_var_id,[i-1,0],[1,length(xGL{i})],uBaseGL{i}');
      netcdf.putVar(ncid,vbase_var_id,[i-1,0],[1,length(xGL{i})],vBaseGL{i}');
      netcdf.putVar(ncid,usurface_var_id,[i-1,0],[1,length(xGL{i})],uSurfaceGL{i}');
      netcdf.putVar(ncid,vsurface_var_id,[i-1,0],[1,length(xGL{i})],vSurfaceGL{i}');
      netcdf.putVar(ncid,umean_var_id,[i-1,0],[1,length(xGL{i})],uMeanGL{i}');
      netcdf.putVar(ncid,vmean_var_id,[i-1,0],[1,length(xGL{i})],vMeanGL{i}');
   end
   netcdf.putVar(ncid,time_var_id,time+dt);

   %Close netcdf
   netcdf.close(ncid)

end

