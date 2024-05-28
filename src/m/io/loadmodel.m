function varargout=loadmodel(path)
%LOADMODEL - load a model using built-in load module 
%
%   check that model prototype has not changed. if so, adapt to new model prototype.
%
%   Usage:
%      md=loadmodel(path)
%      loadmodel path

%check nargout
if nargout>1,
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

try,
	%recover model on file and name it md
	warning off MATLAB:unknownElementsNowStruc;
	warning off MATLAB:load:classNotFound
	struc=load(path,'-mat');
	warning on MATLAB:unknownElementsNowStruc;
	warning on MATLAB:load:classNotFound

	name=char(fieldnames(struc));
	if size(name,1)>1,
		error(['loadmodel error message: file ' path ' contains several variables. Only one model should be present.']); 
	end
	md=struc.(name);
	if nargout,
		varargout{1}=md;
	else
		assignin('caller',name,md);
	end
catch me
	disp(getReport(me))
	error(['could not load model ' path]);
end
