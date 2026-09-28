function vargout = loadmodel(path)
%LOADMODEL - load model from MATLAB file
%
%   Loads the model instance saved in a MATLAB file
%
%   Usage:
%      md=loadmodel(path)
%      loadmodel path

%check nargout
if nargout>1
	error('loadmodel usage error: md=loadmodel(path)');
end

%check existence
if exist(path,'file')
	%do nothing
elseif exist([path '.mat'],'file')
	%add extension
	path = [path '.mat'];
else
	error(['loadmodel error message: file ' path ' does not exist']);
end

try
	%load variables in file 
	warning off MATLAB:unknownElementsNowStruc;
	warning off MATLAB:load:classNotFound
	struc = load(path,'-mat');
	warning on MATLAB:unknownElementsNowStruc;
	warning on MATLAB:load:classNotFound

	%Check that only one variable is present
	name = char(fieldnames(struc));
	if size(name, 1)>1
		error(['loadmodel error message: file ' path ' contains several variables. Only one model should be present.']); 
	end

	%Extract model and make sure it is a model
	md = struc.(name);
	if ~isa(md, 'model')
		error(['variable ''' name ''' saved in ''' path ''' is not of class ''model''']);
	end

	%return model or save in workspace
	if ~nargout
		assignin('caller', name, md);
	else
		vargout(1) = md;
	end

catch me
	disp(getReport(me))
	error(['could not load model from ''' path ''' (see error above)']);
end
