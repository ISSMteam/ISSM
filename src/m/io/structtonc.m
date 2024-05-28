function [var_id,counter] = structtonc(ncid,fieldname,field,depth,var_id,counter,step);
%STRUCTTONC- fill nc file with structure fields
%
%   WARNING: Do not use this function, this function is called
%            by netcdf(model);
%

%update counter
counter   = counter+1;

%Check that field is not empty
if isempty(field) | (isa(field,'struct') & numel(fields(field))==0),
	if(step==1), disp(['skipping ' fieldname ' (empty)...']); end
	return;
end

%Write field class
[var_id,counter] = declareclass(ncid,[fieldname '_class'],class(field),depth,var_id,counter,step);

%Double scalar
if isa(field,'double') & numel(field)==1,
	if step==1,
		var_id(counter) = netcdf.defVar(ncid,fieldname,'NC_DOUBLE',[]);
	else
		netcdf.putVar(ncid,var_id(counter),field);
	end

%Double vector
elseif isa(field,'double') & size(field,2)==1,

if step==1,
		dim_id          = netcdf.defDim(ncid,[fieldname '_size1'],size(field,1));
		var_id(counter) = netcdf.defVar(ncid,fieldname,'NC_DOUBLE',dim_id);
	else
		netcdf.putVar(ncid,var_id(counter),field);
	end

%double matrix
elseif isa(field,'double') & size(field,2)>1,
	if step==1,
		dim1_id         = netcdf.defDim(ncid,[fieldname '_size1'],size(field,1));
		dim2_id         = netcdf.defDim(ncid,[fieldname '_size2'],size(field,2));
		var_id(counter) = netcdf.defVar(ncid,fieldname,'NC_DOUBLE',[dim2_id dim1_id]);
	else
		netcdf.putVar(ncid,var_id(counter),transpose(field));
	end

%string
elseif isa(field,'char') 
		if step==1,
			dim_id          = netcdf.defDim(ncid,[fieldname '_size1'],numel(field));
			var_id(counter) = netcdf.defVar(ncid,fieldname,'NC_CHAR',dim_id);
		else
			netcdf.putVar(ncid,var_id(counter),field);
		end

%Boolean of size 1
elseif isa(field,'logical') & numel(field)==1,
	if step==1,
		var_id(counter) = netcdf.defVar(ncid,fieldname,'NC_BYTE',[]);
	else
		netcdf.putVar(ncid,var_id(counter),int8(field));
	end

%Structures
elseif isa(field,'struct'),
	sublength = numel(field);
	subfields = fields(field);
	allsubfields = '';
	for i=1:length(subfields),
		allsubfields=[allsubfields subfields{i}];
		if i~=length(subfields), allsubfields=[allsubfields ' ']; end
	end
	%Write size
	if step==1,
		var_id(counter) = netcdf.defVar(ncid,[fieldname '_length'],'NC_INT',[]);
	else
		netcdf.putVar(ncid,var_id(counter),sublength);
	end
	counter=counter+1;
	%Write fields
	if step==1,

		dim_id          = netcdf.defDim(ncid,[fieldname '_fields_length'],numel(allsubfields));
		var_id(counter) = netcdf.defVar(ncid,[fieldname '_fields'],'NC_CHAR',dim_id);
	else
		netcdf.putVar(ncid,var_id(counter),allsubfields);
	end
	for n=1:sublength,
		for i=1:length(subfields),
			[var_id,counter] = structtonc(ncid,[fieldname '.' subfields{i} '(' num2str(n) ')'],field(n).(subfields{i}),depth+1,var_id,counter,step);
		end
	end

%Cell
elseif isa(field,'cell'),
	sublength = numel(field);
	%Write size
	if step==1,
		var_id(counter) = netcdf.defVar(ncid,[fieldname '_length'],'NC_INT',[]);
	else
		netcdf.putVar(ncid,var_id(counter),sublength);
	end
	for i=1:sublength,
		[var_id,counter] = structtonc(ncid,[fieldname '{' num2str(i) '}'],field{i},depth+1,var_id,counter,step);
	end

%Objects
elseif isobject(field),
	subfields = fields(field);
	for i=1:length(subfields),
		[var_id,counter] = structtonc(ncid,[fieldname '.' subfields{i}],field.(subfields{i}),depth+1,var_id,counter,step);
	end
else
	disp(['skipping ' fieldname ' (format not supported)...']);
end

function [var_id,counter] = declareclass(ncid,fieldname,field,depth,var_id,counter,step);

if isa(field,'char') 
	if step==1,
		dim_id          = netcdf.defDim(ncid,[fieldname '_length'],numel(field));
		var_id(counter) = netcdf.defVar(ncid,fieldname,'NC_CHAR',dim_id);
	else
		netcdf.putVar(ncid,var_id(counter),field);
	end
else
	error('class name is not a string');
end

%update counter
counter   = counter+1;
