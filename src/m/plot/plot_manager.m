function plot_manager(md,options,subplotwidth,nlines,ncols,i)
%PLOT__MANAGER - distribute the plots, called by plotmodel
%
%   Usage:
%      plot_manager(md,options,subplotwidth,i);
%
%   See also: PLOTMODEL, PLOT_UNIT

%parse options and get a structure of options. 
options=checkplotoptions(md,options);

%get data to be displayed
data=getfieldvalue(options,'data');

%figure out if this is a special plot
if ischar(data),
	switch data,
		case 'boundaries',
			plot_boundaries(md,options,subplotwidth,i);
			return;
		case {'BC','bc'},
			plot_BC(md,options,subplotwidth,i,data);
			return;
		case 'edges'
			plot_edges(md,options,subplotwidth,i,data)
			return
		case 'elementnumbering',
			plot_elementnumbering(md,options,subplotwidth,i);
			return;
		case 'highlightelements',
			plot_highlightelements(md,options,subplotwidth,i);
			return;
		case 'qmumean',
			plot_qmumean(md,options,nlines,ncols,i);
			return;
		case 'qmustddev',
			plot_qmustddev(md,options,nlines,ncols,i);
			return;
		case 'qmuhistnorm',
			plot_qmuhistnorm(md,options,nlines,ncols,i);
			return;
		case 'qmu_mass_flux_segments',
			plot_qmu_mass_flux_segments(md,options,nlines,ncols,i);
			return;
		case 'part_hist',
			plot_parthist(md,options,nlines,ncols,i);
			return;
		case 'part_hist_n',
			plot_parthistn(md,options,nlines,ncols,i);
			return;
		case 'part_hist_w',
			plot_parthistw(md,options,nlines,ncols,i);
			return;
		case 'elements_type',
			plot_elementstype(md,options,subplotwidth,i);
			return;
		case 'vertexnumbering',
			plot_vertexnumbering(md,options,subplotwidth,i);
			return;
		case 'highlightvertices',
			plot_highlightvertices(md,options,subplotwidth,i);
			return;
		case {'basal_drag','basal_dragx','basal_dragy'},
			plot_basaldrag(md,options,subplotwidth,i,data);
			return;
		case 'driving_stress',
			plot_drivingstress(md,options,subplotwidth,i);
			return;
		case 'mesh',
			plot_mesh(md,options,nlines,ncols,i);
			return;
		case 'none',
			if ~exist(options,'overlay'),
				plot_none(md,options,nlines,ncols,i);
				return;
			end
		case 'penalties',
			plot_penalties(md,options,subplotwidth,i);
			return;
		case 'partition',
			plot_partition(md,options,nlines,ncols,i);
			return;
		case 'referential',
			plot_referential(md,options,nlines,ncols,i);
			return;
		case 'riftvel',
			plot_riftvel(md,options,nlines,ncols,i);
			return;
		case 'riftnumbering',
			plot_riftnumbering(md,options,nlines,ncols,i);
			return;
		case 'rifts',
			plot_rifts(md,options,nlines,ncols,i);
			return;
		case 'riftrelvel',
			plot_riftrelvel(md,options,nlines,ncols,i);
			return;
		case 'riftpenetration',
			plot_riftpenetration(md,options,nlines,ncols,i);
			return;
		case 'riftfraction',
			plot_riftfraction(md,options,nlines,ncols,i);
			return;
		case 'sarpwr',
			plot_sarpwr(md,options,subplotwidth,i)
			return
		case 'time_dependant' ,
			plot_vstime(md,options,nlines,ncols,i)
			return
		case 'icefront'
			plot_icefront(md,options,subplotwidth,i,data)
			return
		case 'segments'
			plot_segments(md,options,subplotwidth,i,data)
			return
		case 'quiver'
			data=[md.initialization.vx md.initialization.vy]; %Go ahead and try plot_unit
		case {'strainrate_tensor','strainrate','strainrate_principal','strainrate_principalaxis1','strainrate_principalaxis2','strainrate_principalaxis3',...
				'stress_tensor','stress','stress_principal','stress_principalaxis1','stress_principalaxis2','stress_principalaxis3',...
				'deviatoricstress_tensor','deviatoricstress','deviatoricstress_principal','deviatoricstress_principalaxis1','deviatoricstress_principalaxis2','deviatoricstress_principalaxis3'},
			plot_tensor(md,options,subplotwidth,i,data);
			return;
		case 'thermaltransient_results',
			plot_thermaltransient_results(md,options,subplotwidth,i);
			return;
		case 'transient_movie',
			plot_transient_movie(md,options,subplotwidth,i);
			return;
		case 'transient_results',
			plot_transient_results(md,options,subplotwidth,i);
			return
		case 'transient_field',
			plot_transient_field(md,options,subplotwidth,i);
			return;
	otherwise,
		if isfield(md.results,'TransientSolution') && isfield(md.results.TransientSolution,data)
			plot_transient_movie(md,options,subplotwidth,i);
			return
		elseif ismember(data,properties('model')),
			data=eval(['md.' data ';']);
		else
			error('plot error message: data provided not supported yet. Type plotdoc for help');
		end
	end
end

%Figure out if this is a semi-transparent plot.
if exist(options,'overlay'),
	plot_overlay(md,data,options,nlines,ncols,i);
	return;
end

%Figure out if this is a Google Maps plot.
if exist(options,'googlemaps'),
	plot_googlemaps(md,data,options,nlines,ncols,i);
	return;
end

%Figure out if this is a landsat plot.
if getfieldvalue(options,'landsat',0),
	plot_landsat(md,data,options,nlines,ncols,i);
	return;
end

%Figure out if this is a gridded plot.
if exist(options,'gridded'),
	plot_gridded(md,data,options,nlines,ncols,i);
	return;
end

%Figure out if this is a Section plot
if exist(options,'sectionvalue')
	plot_section(md,data,options,nlines,ncols,i);
	return;
end

%Figure out if this is a Profile plot
if exist(options,'profile')
	plot_profile(md,data,options,nlines,ncols,i);
	return;
end

%process data and model
[x y z elements is2d isplanet]=processmesh(md,data,options);
[data2 datatype]=processdata(md,data,options);

%standard plot:
if exist(options,'asymsubplot')
	id=getfieldvalue(options,'asymsubplot',i);
	subplotmodel(nlines,ncols,id,options);
else
	subplotmodel(nlines,ncols,i,options);
end

% clear subplot option
if exist(options,'subplot')
    cla;
end

%plot unit
plot_unit(x,y,z,elements,data2,is2d,isplanet,datatype,options);

%apply all options
if datatype==3,
	options=changefieldvalue(options,'colorbar',2);
	if exist(options,'contourlevels'),
		data2=data;
	end
end

applyoptions(md,data2,options);

%do ground overlay on kml plot_unit? 
if (strcmpi(getfieldvalue(options,'kmlgroundoverlay','off'),'on')),
	if ((nlines*ncols~=1) | (i~=1)),
		error('cannot kmlgroundoverlay on multi-plots');
	end

	%call routine to build kml file and image that goes with it.
	kmlgroundoverlay(md,options);
end
