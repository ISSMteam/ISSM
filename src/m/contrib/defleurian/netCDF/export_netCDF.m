function export_netCDF(md,filename)
%verbosity of the code, 0 is no messages, 5 is chatty
	verbose = 5;
	if exist(filename),
		delete(filename)
		% disp(sprintf('File %s allready exist', filename));
		% prompt = 'Give a new name or "delete" to replace: ';
		% newname = input(prompt,'s');
		% if strcmp(newname,'delete')
		% 	delete(filename)
		% else
		% 	disp(sprintf('New file name is %s ', newname));
		% 	filename=newname
		% end
	end
	%open file and write description
	mode = netcdf.getConstant('NC_NETCDF4');
	mode = bitor(mode,netcdf.getConstant('NC_NOCLOBBER')); %NOCLOBBER to avoid overwrite
	ncid = netcdf.create(filename,mode);
	netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'description',['Results for run ' md.miscellaneous.name]);
	netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'history',['Created ' datestr(now)]);

	%gather geometry and timestepping as dimensions
	if isempty(fieldnames(md.results)),
		%results as no field so no time is present
		Duration = 0;
	else
		resfields = fieldnames(md.results);
		Duration = size(eval(['md.results. ' resfields{1} ]),2);
	end
	if Duration>0,
		StepNum = Duration;
	else
		StepNum=1;
	end
	DimSize(1).index=netcdf.defDim(ncid,'Time',StepNum);  %time is the first dimension
	[DimSize(1).name,DimSize(1).value]=netcdf.inqDim(ncid,DimSize(1).index);
	DimValue(1)=DimSize(1).value;
	DimSize(2).index=netcdf.defDim(ncid,'UnLim',netcdf.getConstant('NC_UNLIMITED')); % we add an unlimited dimension if needed
	[DimSize(2).name,DimSize(2).value]=netcdf.inqDim(ncid,DimSize(2).index);
	DimValue(2)=DimSize(2).value;
	% adding mesh related dimensions
	dimlist=[2 40 md.mesh.numberofelements md.mesh.numberofvertices size(md.mesh.elements,2)];
	dimnames=["DictDummy" "StringLength" "EltNum" "VertNum" "VertPerElt"];
	if isprop(md.mesh, 'edges'),
		dimlist(end+1)=md.mesh.numberofedges;
		dimnames(end+1)="EdgeNum";
	else
		dimlist(end+1)=0;
		dimnames(end+1)="EdgeNum";
	end
	if verbose > 0,
		disp('===Creating dimensions ===');
	end
	%define netcdf dimensions
	for i=1:length(dimlist)
		% do not add the dimension if it exists already
		if sum(dimlist(i) == DimValue) == 0
			DimSize(i+2).index=netcdf.defDim(ncid,dimnames(i),dimlist(i));
			[DimSize(i+2).name,DimSize(i+2).value]=netcdf.inqDim(ncid,DimSize(i+2).index);
			DimValue(i+2)=DimSize(i+2).value;
		end
	end
	issmclasses = fieldnames(md)';
	typelist={'half', 'single','double','int8','int16'...
		  ,'int32','int64','uint8','uint16','uint32'...
		  ,'uint64','logical','char','string'};  %all malab types that are 0D

	for cl=1:length(issmclasses),
		if isempty(md.(issmclasses{cl})),
			disp(sprintf("md.%s is empty and will be left as default",issmclasses{cl}));
		else
			subclasses=fieldnames(md.(issmclasses{cl}))';
			for sc=1:length(subclasses),
				if sum(strcmp(class(md.(issmclasses{cl}).(subclasses{sc})), typelist)) == 0,
					issmclasses = [issmclasses class(md.(issmclasses{cl}).(subclasses{sc}))];
				end
			end
		end
	end
	%get all model classes and create respective groups
	groups=fieldnames(md);
	if verbose > 0,
		disp('===Creating and populating groups===');
	end
	for i=1:length(groups),
		if verbose >1,
			disp(sprintf('===Now treating %s===',groups{i}));
		end
		if any(strcmp(groups{i}, {'qmu'})),
			disp('qmu is skipped until it is more stable');
			continue
		end
		if any(strcmp(groups{i},{'radaroverlay'})),
			disp(sprintf('%s is skipped.',groups{i}));
			continue
		end
		groupID=netcdf.defGrp(ncid,groups{i});
		%In each group gather the fields of the class
		if isempty(md.(groups{i})),
			disp(sprintf("WARNING: md.%s is empty, we skip it.",groups{i}))
			continue
		end
		fields=fieldnames(md.(groups{i}));
		if isempty(fields),
			disp(sprintf("WARNING: md.%s as no fields, we skip it.",groups{i}))
			continue
		end
		%looping on fields in each group
		for j=1:length(fields),
			Var=md.(groups{i}).(fields{j});
			%treatment for lists
			if isa(Var,'cell')
				Stdlist=false;  %first assume it is not a standard list
				if length(Var) == 0
					Stdlist=true;  %It is empty and so standard (for us)
				else
					for k=1:length(typelist)
						if isa(Var{1},typelist{k})
							Stdlist=true;  %if the list is of a known type (to matlab) if not it is probably some exotic ISSM stuff
						end
					end
				end
				%print the issm class as a classtype attribute
				klass = class(md.(groups{i}));
				klasstring = strcat(klass, '.',klass);
				netcdf.putAtt(groupID,netcdf.getConstant('NC_GLOBAL'),'classtype',klasstring);
				if(Stdlist)  % this is a standard or empty list just proceed
					if verbose > 4,
						disp(sprintf("=££=creating var for %s.%s with classtype : %s",groups{i}, fields{j}, klasstring))
					end
					if ~isempty(Var) && isa(Var{1}, 'char'),  % we have a char array, pad it to a given length
						Var=char(Var)';
					end
					[DimSize,DimValue,varid]=CreateVar(ncid,Var,groupID,fields{j},DimSize,DimValue);
					if ~isempty(varid),
						FillVar(Var,groupID,varid);
					end

				else % this is a list of fields, specific treatment needed (perhaps)
					if verbose > 4,
						disp(sprintf("=??=we have a list of fields for %s.%s with classtype : %s",groups{i}, fields{j}, klasstring));
					end
					if strcmp(groups{i}, 'outputdefinition'),
						listsize=length(Var);
						for k=1:listsize,
							subgroupname=md.(groups{i}).(fields{j}){k}.definitionstring;
							subgroupID=netcdf.defGrp(groupID,subgroupname);
							klass=class(md.(groups{i}).(fields{j}){k});
							klasstring = strcat(klass, '.',klass);
							netcdf.putAtt(subgroupID,netcdf.getConstant('NC_GLOBAL'),'classtype',klasstring);
							subfields=fieldnames(md.(groups{i}).(fields{j}){k});
							for l=1:length(subfields)
								if verbose > 4,
									disp(sprintf("=--=creating var for %s.%s[%i].%s",groups{i}, fields{j}, k, subfields{l}));
								end
								Var = md.(groups{i}).(fields{j}){k}.(subfields{l});
								if sum(numel(Var) == size(Var)) == 0,  %this is a 2D array or more (and not a vector with dimension 2 = 1)
									Var = Var';
								end
								[DimSize,DimValue,varid]=CreateVar(ncid,Var,subgroupID,subfields{l},DimSize,DimValue);
								if ~isempty(varid),
									FillVar(Var,subgroupID,varid);
								end
							end
						end
					else
						disp(sprintf("WARNING: unknown treatment for md.%s",groups{i}));
					end
				end
			elseif sum(strcmp(class(Var), typelist))==1, %this is a standard matlab class with no subgrouping
				if verbose > 4,
					disp(sprintf("====creating var for %s.%s", groups{i}, fields{j}))
				end
				klass=class(md.(groups{i}));
				klasstring = strcat(klass, '.',klass);
				netcdf.putAtt(groupID,netcdf.getConstant('NC_GLOBAL'),'classgroup',groups{i});
				netcdf.putAtt(groupID,netcdf.getConstant('NC_GLOBAL'),'classtype',klasstring);
				if sum(numel(Var) == size(Var)) == 0,  %this is a 2D array or more (and not a vector with dimension 2 = 1)
					Var = Var';
				end

				[DimSize,DimValue,varid]=CreateVar(ncid,Var,groupID,fields{j},DimSize,DimValue);
				if ~isempty(varid),
					FillVar(Var,groupID,varid);
				end


			elseif isa(Var,'struct')  % structures need special treatment
				if strcmp(groups{i}, 'results'),
					klasstring=strcat(groups{i} ,'.', groups{i});
					netcdf.putAtt(groupID,netcdf.getConstant('NC_GLOBAL'),'classtype',klasstring);
					Listsize= length(md.(groups{i}).(fields{j}));
					subgroupname=fields{j};
					subgroupID=netcdf.defGrp(groupID,subgroupname);
					klasstring='results.solutionstep';
					netcdf.putAtt(subgroupID,netcdf.getConstant('NC_GLOBAL'),'classtype',klasstring);
					subfields=fieldnames(md.(groups{i}).(fields{j}));
					if isempty(subfields),
						disp(sprintf("WARNING: md.%s.%s as no subfields, we skip it.",groups{i}, fields{j}));
						continue
					end
					for k=1:length(subfields),
						if ~ismember(subfields{k}, {'errlog', 'outlog', 'SolutionType'})
							StackedVar=restable();
							for l=1:Listsize,
								Var = md.(groups{i}).(fields{j})(l).(subfields{k});
								if length(Var) == 0,
									%Some variables only have data on the first step
									break
								end
								lastindex=l;
								StackedVar=StackedVar.update(Var);
							end
							if verbose > 4,
								disp(sprintf("=$$=creating var for %s.%s.%s",groups{i}, fields{j}, subfields{k}));
								disp(sprintf("last index on the list is %i",lastindex));
							end
							StackedVar=StackedVar.finalize(lastindex);
							[DimSize,DimValue,varid]=CreateVar(ncid,StackedVar,subgroupID,subfields{k},DimSize,DimValue);
							if ~isempty(varid),
								FillVar(StackedVar,subgroupID,varid);
							end
						elseif ismember(subfields{k}, {'SolutionType'})
							%We just add solution type once as an attribute
							Var = md.(groups{i}).(fields{j})(1).(subfields{k});
							[DimSize,DimValue,varid]=CreateVar(ncid,Var,subgroupID,subfields{k},DimSize,DimValue);
							if ~isempty(varid),
								FillVar(Var,subgroupID,varid);
							end

						end
					end
				elseif strcmp(groups{i}, 'toolkits'),
					klasstring=strcat(groups{i} ,'.', groups{i});
					netcdf.putAtt(groupID,netcdf.getConstant('NC_GLOBAL'),'classtype',klasstring);
					if verbose > 4,
						disp(sprintf("=}{=creating var for %s.%s",groups{i}, fields{j}));
					end

					[DimSize,DimValue,varid]=CreateVar(ncid,Var,groupID,fields{j},DimSize,DimValue);
					if ~isempty(varid),
						FillVar(Var,groupID,varid);
					end

				elseif isempty(fieldnames(md.(groups{i}).(fields{j}))) % this is an empty struct, jus treat it as normal
					klass=class(md.(groups{i}));
					klasstring = strcat(klass, '.',klass);
					netcdf.putAtt(groupID,netcdf.getConstant('NC_GLOBAL'),'classtype',klasstring);
					if verbose > 4,
						disp(sprintf("=[]=creating var for %s.%s",groups{i}, fields{j}));
					end

					[DimSize,DimValue,varid]=CreateVar(ncid,Var,groupID,fields{j},DimSize,DimValue);
					if ~isempty(varid),
						FillVar(Var,groupID,varid);
					end

				else
					disp(sprintf("WARNING, md.%s.%s is not treated as it does not fall in one of the existing cases with class '%s'.",groups{i}, fields{j}, class(md.(groups{i}).(fields{j}))))
				end
			elseif sum(strcmp(class(Var), issmclasses)) == 1,  % that is an issm class
				if strcmp(class(Var), 'solution'),
					if verbose > 4,
						disp(sprintf("=$$=creating var for %s.%s",groups{i}, fields{j}))
						disp("NEED treatment")
					end
				elseif strcmp(class(Var), 'dict'),  %we have potential for a dict in py not to sure what it translates to here.
					if verbose > 4,
						disp(sprintf("=WW=creating var for %s.%s",groups{i}, fields{j}))
						disp("NEED Treatment")
					end

				else
					klass=class(md.(groups{i}));
					klasstring = strcat(klass, '.',klass);
					netcdf.putAtt(groupID,netcdf.getConstant('NC_GLOBAL'),'classtype',klasstring);
					subgroupID=netcdf.defGrp(groupID,fields{j});
					klass=class(md.(groups{i}).(fields{j}));
					klasstring = strcat(klass, '.',klass);
					netcdf.putAtt(subgroupID,netcdf.getConstant('NC_GLOBAL'),'classtype',klasstring);
					subfields=fieldnames(Var);
					for k=1:length(subfields),
						if sum(strcmp(subfields{k},["outlog" "errlog"])) == 0,
							if verbose > 4,
								disp(sprintf("+==+creating var for %s.%s.%s",groups{i}, fields{j}, subfields{k}))
							end
							Var=md.(groups{i}).(fields{j}).(subfields{k});
							[DimSize,DimValue,varid]=CreateVar(ncid,Var,subgroupID,subfields{k},DimSize,DimValue);
							if ~isempty(varid),
								FillVar(Var,subgroupID,varid);
							end
						end
					end

				end
			else
				disp(sprintf("WARNING, md.%s.%s is not treated as it does not fall in one of the existing cases with class '%s'.",groups{i}, fields{j}, class(Var)))
			end
		end
	end
	netcdf.close(ncid);
end

function [DimSize,DimValue,varid]=CreateVar(ncid,Var,groupID,field,DimSize,DimValue)
% Grab dimensions
	varsize=size(Var);
	varlength=length(Var);
	% treating scalar string or bool as atribute
	if isa(Var,'logical'),
		if Var,
			LogicString='True';
		else,
			LogicString='False';
		end
		netcdf.putAtt(groupID,netcdf.getConstant('NC_GLOBAL'),field,LogicString);
		varid=[];

	elseif isa(Var,'char'),
		if strcmp(field,'name'),  % it looks like netCDF does not like attributes that are called "name"
			field = 'varname';
		end
		if size(Var,1) <= 1  %that is a single string or empty
			netcdf.putAtt(groupID,netcdf.getConstant('NC_GLOBAL'),field,Var);
			varid=[];
		else  % that is a character array
			[dims,DimSize,DimValue]=GetDims(ncid,Var,DimSize,DimValue);
			varid = netcdf.defVar(groupID,field,'NC_CHAR',dims);
			if numel(Var)>1
				netcdf.defVarDeflate(groupID,varid,true,true,4);
			end
		end

	elseif isa(Var,'double'), %dealing with arrays
		if all(mod(Var, 1) == 0, 'all')  %those are actually integers,
			[dims,DimSize,DimValue]=GetDims(ncid,Var,DimSize,DimValue);
			varid = netcdf.defVar(groupID,field,'NC_INT64',dims);
			if numel(Var)>1
				netcdf.defVarDeflate(groupID,varid,true,true,4);
			end
		else
			[dims,DimSize,DimValue]=GetDims(ncid,Var,DimSize,DimValue);
			varid = netcdf.defVar(groupID,field,'NC_DOUBLE',dims);
			if numel(Var)>1
				netcdf.defVarDeflate(groupID,varid,true,true,4);
			end
		end
	elseif isa(Var,'cell'),
		% cells can be a range of things, what are we dealing with here
		if isempty(Var),
			netcdf.putAtt(groupID,netcdf.getConstant('NC_GLOBAL'),field,'emptycell');
			varid=[];
		else
			[dims,DimSize,DimValue]=GetDims(ncid,Var,DimSize,DimValue);
			if isa(Var{1}, 'double'),
				varid = netcdf.defVar(groupID,field,'NC_DOUBLE',dims);
				if numel(Var)>1
					netcdf.defVarDeflate(groupID,varid,true,true,4);
				end
			else
				varid = netcdf.defVar(groupID,field,'NC_CHAR',dims);
				if numel(Var)>1
					netcdf.defVarDeflate(groupID,varid,true,true,4);
				end
			end
		end
	elseif isa(Var,'struct'),
		if isempty(fieldnames(Var)),
			netcdf.putAtt(groupID,netcdf.getConstant('NC_GLOBAL'),field,'emptystruct');
			varid=[];
		else
			%Start by getting the structure fields and size
			[dims,DimSize,DimValue]=GetDims(ncid,Var,DimSize,DimValue);
			varid = netcdf.defVar(groupID,field,'NC_CHAR',dims);
			if numel(Var)>1
				netcdf.defVarDeflate(groupID,varid,true,true,4);
			end
		end
	else
		disp(sprintf('no support for class %s of field %s',class(Var),field));
		varid=[];
	end
	return
end


function FillVar(Var,groupID,varid)
% Grab dimensions
	varsize=size(Var);
	varlength=length(Var);
	% treating scalar string or bool as atribute
	if isa(Var,'double'), %dealing with arrays
		if all(mod(Var, 1) == 0, 'all')  %those are actually integers,
			Var = int64(Var);
		end
		if length(Var)==0,
			netcdf.putVar(groupID,varid,NaN);
		else
			netcdf.putVar(groupID,varid,Var);
		end
	elseif isa(Var,'char'),  % at this point this should be a character array
		netcdf.putVar(groupID,varid,Var);
	elseif isa(Var,'cell'),  % there can be a number of things in a cell array
		for i=1:length(Var),
			if isa(Var{i},'char')  %for characters we limit the size to 40 for now
				if length(Var)>1,
					count=[min(length(Var{i}),40), 1];
					startpoint=[0 i-1];
				else
					count=min(length(Var{i}),40);
					startpoint=0;
				end

				if length(Var{i})>40,
					netcdf.putVar(groupID,varid,startpoint,count,Var{i}(1:40));
					disp(sprintf('some variable have been truncated'));
				else
					netcdf.putVar(groupID,varid,startpoint,count,Var{i});
				end
			elseif isa(Var{i},'double')
				startpoint=[i-1];
				count=[1 length(Var{i}) ndims(Var{i})];
				for j=1:ndims(Var{i}),
					startpoint=[startpoint 0];
				end
				netcdf.putVar(groupID,varid,startpoint,count,Var{i});
			else
				disp(sprintf("WARNING: cell of class %s is not supported.",class(Var{i})))
			end
		end
	elseif isa(Var,'struct'),
		%Start by getting the structure fields and size
		locfields=fieldnames(Var);
		for i=1:length(locfields),
			for j=1:2,
				if j==1,
					CharVar=locfields{i}';
					disp(size(CharVar))
					if length(CharVar)==0
						CharVar='emptystruct';
					end
					startpoint=[0,0,i-1]
				else
					if isa(Var.(locfields{i}),'char'),
						CharVar=Var.(locfields{i})';
					else
						CharVar=num2str(Var.(locfields{i}))';
					end
					if length(CharVar)==0
						CharVar='emptystruct';
					end
					startpoint=[0,1,i-1]
				end

				extent=[min(length(CharVar),40), 1, 1]
				if length(CharVar)>40,
					netcdf.putVar(groupID,varid,startpoint,extent,CharVar(1:40));
					disp(sprintf('some variable have been truncated'));
				else
					netcdf.putVar(groupID,varid,startpoint,extent,CharVar);
				end
			end
		end
	else
		disp(sprintf('no support for class %s',class(Var)));
	end
	return
end

function [dims,DimSize,DimValue]=GetDims(ncid,Var,DimSize,DimValue)
	dims=[];
	celldims=[];
	dim=ndims(Var);
	if isa(Var,'struct'),
		varsize=length(fieldnames(Var));
	else
		varsize=size(Var);
		if isa(Var, 'cell')
			%we add the dimension of the cells themselves,
			%that will most probably fail if cells have different sizes
			for i=1:dim,
				newdim=size(Var{i});
				if ~ismember(newdim, celldims),
					celldims=[celldims newdim];
				end
			end
		end
	end
	varsize=[varsize celldims];
	alldim=length(varsize);
	if dim>0,
		for i=1:alldim,
			if size(Var, i)>1 || i>dim || isa(Var, 'struct'),  %we skip dimensions with zero lenght but want to add dimensions from cells
				indsize=find(varsize(i)==DimValue);
				if length(indsize)>0
					dims=[dims DimSize(indsize).index];
				else
					indsize=length(DimSize)+1;
					DimSize(indsize).index=netcdf.defDim(ncid,['DimNum' num2str(indsize)],varsize(i));
					[DimSize(indsize).name,DimSize(indsize).value]=netcdf.inqDim(ncid,DimSize(indsize).index);
					DimValue(indsize)=DimSize(indsize).value;
					dims=[dims DimSize(indsize).index];
				end
			end
		end
	end
	if isa(Var, 'cell') && isa(Var{1}, 'char'),
		%if we have an cell variable with strings we need to add a stringlength
		dims=[dims DimSize(4).index];
	end
	% struct also need an extra dimension 2, but only if non empty
	if isa(Var,'struct'),
		dims=[DimSize(4).index DimSize(3).index, dims];
	end
end
