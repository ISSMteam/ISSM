function ofname = interpXarrayGridToMesh(fname, X, Y, varargin)
	%{
	% Explain
	%
	% Usage
	% -----
	% md = loadmodel('Models/Mesh.mat');
	% data = interpXarrayGridToMesh(nc,md.mesh.x,md.mesh.y,'xname','x','yname','y','varname','variable');
	%
	% Inputs
	%  fname  - input file
	%
	% Outputs
	% -------
	%  ofname - return nc file.
	%}

	% Get options
	options   = pairoptions(varargin{:});
	isverbose = getfieldvalue(options,'verbose',0);
	ismethod  = getfieldvalue(options,'method','linear');
	varname   = getfieldvalue(options,'varname','');
	xname   = getfieldvalue(options,'xname','X');
	yname   = getfieldvalue(options,'yname','Y');

	% Initialize temporal file.
	if ~exist('./temp')
		mkdir temp
	end
	inputgrid = [tempname('./temp'), '.mat'];
	ofname = [tempname('./temp/'), '.nc'];

	if isverbose
		fprintf('   Inputgrid: %s\n', inputgrid);
		fprintf('   Ouputname: %s\n', ofname);
	end

	% Initialize function directory
	[func_dirname, ~, ~] = fileparts(which('interpXarrayGridToMesh'));

	% Prepare the xy grid.
	if 0
		fid = fopen(inputgrid, 'w');
		for i = 1:length(X)
			fprintf(fid, '%f %f\n', X(i), Y(i));
		end
		fclose(fid);
	else
		disp('   Do interpolation');
		save(inputgrid,'-v7.3','X','Y');
	end

	% Do interpolation
	command = ['python3 ' func_dirname '/interpXarrayGridToMesh.py -xname ' xname ' -yname ' yname];
	if isverbose
		command = [command ' -v 1'];
	end
	command = [command ' -inputformat mat'];
	command = [command ' -method ' ismethod ' ' fname ' ' varname ' ' inputgrid ' ' ofname];

	system(command);

	if ~exist(ofname)
		disp(command);
		error(sprintf('ERROR: We cannot find %s', ofname));
	end
end


