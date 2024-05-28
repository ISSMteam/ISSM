function []=convertmattoarch(matfile,archfile);
%    convertmattoarch -- function to convert mat-format file to arch-format file
%
%    usage:
%        convertmattoarch(matfile,archfile);
%    where:
%        matfile		name of mat-format file
%        archfile		name of arch-format file (optional)
%

	if ~exist('matfile','var') || isempty(matfile)
		help convertmattoarch
		error('convertmattoarch usage error.');
	end

	[pathstr,name,ext]=fileparts(matfile);
	if isempty(ext)
		ext='.mat';
	end
	matfile=fullfile(pathstr,[name ext]);

	if ~exist('archfile','var') || isempty(archfile)
		archfile=fullfile(pathstr,[name '.arch']);
	end

	if exist(archfile,'file')
		delete(archfile);
	end

	a=load(matfile,'-mat');
	disp(sprintf('mat-format file ''%s'' read.',matfile));
	fnames=fieldnames(a);

	for i=1:length(fnames)
		if isstruct(a.(fnames{i})) || iscell(a.(fnames{i}))
			warning('field ''%s'' is of class ''%s'' and will not be written.',fnames{i},class(a.(fnames{i})));
		else
			% matlab writes the dimensions reversed and matrices transposed into binary files, so compensate for that
			archwrite(archfile,fnames{i},transpose(a.(fnames{i})));
			disp(sprintf('field ''%s'' of class ''%s'' and size [%dx%d] written.',...
			             fnames{i},class(a.(fnames{i})),size(a.(fnames{i}),1),size(a.(fnames{i}),2)));
		end
	end
	disp(sprintf('arch-format  file ''%s'' written.',archfile));

end

