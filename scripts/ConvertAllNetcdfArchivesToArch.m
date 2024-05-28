root=pwd;
cd('../test/Archives/');
ncdir=pwd;
cd(root);

file_regexp=strcat(ncdir,'/*.nc');
flist=dir(file_regexp);
fnames={flist.name}';

pattern='.nc';
replacement='.arch';
for i=1:numel(fnames)
	% Remove the extension '.nc' and replace with '.arch'
	arch_file_name=regexprep(fnames{i},pattern,replacement);
	disp(['Converting ' fnames{i} ' to .arch format']);
	% Reconstruct file path
	fpath=strcat(ncdir,'/',fnames{i});
	archpath=strcat(ncdir,'/',arch_file_name);
	netcdf2arch(fpath,archpath);
end
