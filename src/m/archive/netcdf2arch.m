function netcdf2arch(netcdf_filename,arch_filename) % {{{
%NETCDF2ARCH - Convert a netcdf file into an archive file.
%
%	Usage:
%		netcdf2arch('Archive101.nc','Archive101.arch');
	ncdata=readnetcdf(netcdf_filename);
	ncstruct=struct();

	variables = ncdata.VarArray;
	for i=1:size(variables,2),
		fieldname=deblank(variables(i).Str);
		fieldvalue=double(squeeze(variables(i).Data));
		ncstruct.(fieldname)=fieldvalue;
	end

	variables=ncdata.AttArray;
	for i=1:size(variables,2),
		fieldname=deblank(variables(i).Str);
		fieldvalue=double(variables(i).Val);
		ncstruct.(fieldname)=fieldvalue;
	end
	ncfields=fieldnames(ncstruct);
	% first, remove old arch file
	delete(arch_filename);
	% write data to file 
	for i=1:numel(ncfields)
		archwrite(arch_filename,ncfields{i},ncstruct.(ncfields{i}));
	end
end%}}}
