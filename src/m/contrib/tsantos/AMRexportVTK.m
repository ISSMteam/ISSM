function AMRexportVTK(filename,model,varargin)
%AMREXPORTVTK -  vtk export
%
%   function AMRexportVTK(filename,model)
%   creates a directory with the vtk files for displays in paraview
%   (only work for triangle based on their number of nodes)
%   By default only the results are exported
%
%   USAGE:
%      AMRexportVTK(filename,model,varargin)
%
%   EXAMPLE:
%      AMRexportVTK('ResultSimulation1',md)

[path,name,ext]=fileparts(filename);
separator=filesep;
mkdir(filename);

%get the element related variables
if dimension(model.mesh)==2,
	points=[model.mesh.x model.mesh.y zeros(model.mesh.numberofvertices,1)];
else
	error('Dimension not supported yet \n');
end

[num_of_points,dim]=size(points);
[num_of_elt]=size(model.mesh.elements,1);
[point_per_elt]=size(model.mesh.elements,2);

%Select the type of element function of the number of nodes per elements
if point_per_elt==3;
	celltype=5; %triangles
elseif point_per_elt==6;
	error('Element type not supported yet \n')
else
	error('Your Element definition is not taken into account \n');
end

%this is the result structure (just the Transient solution)
res_struct=model.results;
%checking for results
if (length(fields(res_struct))>0);
	%Getting all the solutions of the model
	solnames=fields(res_struct);
	num_of_sols=length(solnames);
	num_of_timesteps=1;
	%building solution structure 
	for i=1:num_of_sols
		sol_struct{i}=res_struct.(solnames{i});
		%looking for multiple time steps
		if(size(sol_struct{i},2)>num_of_timesteps);
			num_of_timesteps=size(sol_struct{i},2);
			outstep=model.timestepping.time_step*model.settings.output_frequency;
	  end
  end
else
	num_of_timesteps=1;
end
for step=1:num_of_timesteps;
	
	timestep=step;
 	if(isfield(model.results.TransientSolution,'MeshElements'))
   	index = model.results.TransientSolution(step).MeshElements;
   	x     = model.results.TransientSolution(step).MeshX;
   	y     = model.results.TransientSolution(step).MeshY;
   else
   	index = model.mesh.elements;
   	x     = model.mesh.x;
   	y     = model.mesh.y;
   end
	
	points=[x y zeros(size(x))];
	[num_of_points,dim]=size(points);
	[num_of_elt]=size(index,1);

	fid = fopen(strcat(path,filesep,name,filesep,'timestep.vtk',int2str(timestep),'.vtk'),'w+');
	fprintf(fid,'# vtk DataFile Version 2.0 \n');
	fprintf(fid,'Data for run %s \n',model.miscellaneous.name);
	fprintf(fid,'ASCII \n');
	fprintf(fid,'DATASET UNSTRUCTURED_GRID \n');
	
	fprintf(fid,'POINTS %d float\n',num_of_points);
	if(dim==3);
		s='%f %f %f \n';
	elseif(dim==2);
		s='%f %f \n';
  end
	P=[points zeros(num_of_points,3-dim)];
	fprintf(fid,s,P');
	
	fprintf(fid,'CELLS %d %d\n',num_of_elt,num_of_elt*(point_per_elt+1));
	s='%d';
	for j=1:point_per_elt
		s=horzcat(s,{' %d'});
  end
	s=cell2mat(horzcat(s,{'\n'}));
		fprintf(fid,s,[(point_per_elt)*ones(num_of_elt,1) index-1]');
	
	fprintf(fid,'CELL_TYPES %d\n',num_of_elt);
	s='%d\n';
	fprintf(fid,s,celltype*ones(num_of_elt,1));
	fprintf(fid,'POINT_DATA %s \n',num2str(num_of_points));

	%loop over the different solution structures
	if (exist('num_of_sols'));
		for j=1:num_of_sols
			%dealing with results on different timesteps
			if(size(sol_struct{j},2)>timestep);
				timestep = step;
			else
				timestep = size(sol_struct{j},2);
			end
			
			%getting the number of fields in the solution
			fieldnames=fields(sol_struct{j}(timestep));
			num_of_fields=length(fieldnames);
			
			%check which field is a real result and print
			for k=1:num_of_fields
				if ((numel(sol_struct{j}(timestep).(fieldnames{k})))==num_of_points);
					%paraview does not like NaN, replacing
					nanval=find(isnan(sol_struct{j}(timestep).(fieldnames{k})));
					sol_struct{j}(timestep).(fieldnames{k})(nanval)=-9999;
					%also checking for verry small value that mess up
					smallval=(abs(sol_struct{j}(timestep).(fieldnames{k}))<1.0e-20);
					sol_struct{j}(timestep).(fieldnames{k})(smallval)=0.0;
					fprintf(fid,'SCALARS %s float 1 \n',fieldnames{k});
					fprintf(fid,'LOOKUP_TABLE default\n');
					s='%e\n';
					fprintf(fid,s,sol_struct{j}(timestep).(fieldnames{k}));
				end
			end
			fprintf(fid,'CELL_DATA %s \n',num2str(num_of_elt));
			for k=1:num_of_fields
				if ((numel(sol_struct{j}(timestep).(fieldnames{k})))==num_of_elt);
					%paraview does not like NaN, replacing
					nanval=find(isnan(sol_struct{j}(timestep).(fieldnames{k})));
					sol_struct{j}(timestep).(fieldnames{k})(nanval)=-9999;
					%also checking for verry small value that mess up
					smallval=(abs(sol_struct{j}(timestep).(fieldnames{k}))<1.0e-20);
					sol_struct{j}(timestep).(fieldnames{k})(smallval)=0.0;
					fprintf(fid,'SCALARS %s float 1 \n',fieldnames{k});
					fprintf(fid,'LOOKUP_TABLE default\n');
					s='%e\n';
					fprintf(fid,s,sol_struct{j}(timestep).(fieldnames{k}));
				end		
			end
	  end
  end
	%loop on arguments, if something other than result is asked, do
	%it now
	for j= 1:nargin-2
		res_struct=model.(varargin{j});
		fieldnames=fields(res_struct);
		num_of_fields=length(fieldnames);
		for k=1:num_of_fields
			if ((numel(res_struct.(fieldnames{k})))==num_of_points);
				%paraview does not like NaN, replacing
				nanval=find(isnan(res_struct.(fieldnames{k})));
				res_struct.(fieldnames{k})(nanval)=-9999;
				%also checking for verry small value that mess up
				smallval=(abs(res_struct.(fieldnames{k}))<1.0e-20);
				res_struct.(fieldnames{k})(smallval)=0.0;
				fprintf(fid,'SCALARS %s float 1 \n',fieldnames{k});
				fprintf(fid,'LOOKUP_TABLE default\n');
				s='%e\n';
				fprintf(fid,s,res_struct.(fieldnames{k}));
				%check for forcings	
			elseif (size(res_struct.(fieldnames{k}),1)==num_of_points+1);
				%paraview does not like NaN, replacing
				nanval=find(isnan(res_struct.(fieldnames{k})));
				res_struct.(fieldnames{k})(nanval)=-9999;
				%also checking for verry small value that mess up
				smallval=(abs(res_struct.(fieldnames{k}))<1.0e-20);
				res_struct.(fieldnames{k})(smallval)=0.0;
				if (size(res_struct.(fieldnames{k}),2)==num_of_timesteps),
					fprintf(fid,'SCALARS %s float 1 \n',fieldnames{k});
					fprintf(fid,'LOOKUP_TABLE default\n');
					s='%e\n';
					fprintf(fid,s,res_struct.(fieldnames{k})(1:end-1,timestep));
				else,
					%forcing and results not on the same timestep,need some treatment
					fprintf(fid,'SCALARS %s float 1 \n',fieldnames{k});
					fprintf(fid,'LOOKUP_TABLE default\n');
					index=1;
					currenttime=((timestep-1)*outstep)+model.timestepping.start_time;
					while (res_struct.(fieldnames{k})(end,index)<=currenttime);
						if index==size(res_struct.(fieldnames{k}),2)
							break
						end	
						index=index+1;
		      end
					uptime=res_struct.(fieldnames{k})(end,index);
					uplim=res_struct.(fieldnames{k})(1:end-1,index);
					while (res_struct.(fieldnames{k})(end,index)>=currenttime);
						if index==1
							break
			      end
						index=index-1;
		      end
					lowtime=res_struct.(fieldnames{k})(end,index);
					lowlim=res_struct.(fieldnames{k})(1:end-1,index);
					if uptime==currenttime,
						interp=uplim;
					elseif lowtime==currenttime,
						interp=lowlim;
					else
						interp=lowlim+(uplim-lowlim)*((currenttime-lowtime)/(uptime-lowtime));
					end
					s='%e\n';
					fprintf(fid,s,interp);
				end
		  end		
		end 
	end
	fclose(fid);
end
