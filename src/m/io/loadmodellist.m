function varargout=loadmodellist(path)
%LOADMODELLIST - load a list of models from a .mat file
%
%   Loads the modellist instance saved in a .mat file, converting from an
%   older structure format to a modellist object if necessary.
%
%   Usage:
%      mds=loadmodellist(path)
%      loadmodellist path

%check nargout
if nargout>1
	error('loadmodellist usage error: mds=loadmodellist(path)');
end
%check existence
if ~exist(path)
	error(['loadmodellist error message: file ' path ' does not exist']);
end

%check that the file is readable
[stat,mess]=fileattrib(path);
if( stat==0 | mess.UserRead~=1),
	error(['loadmodellist error message: file ' path ' is not readable (permission dinied).']);
end

%check number of variables
if length(whos('-file',path))>1
	error(['loadmodellist error message: file ' path ' contains several variables. Only one model should be present.']);
end

try,
	struc=load(path,'-mat');

	%get name of model variable
	fieldname=char(fieldnames(struc));
	mds=eval(['struc.' fieldname]);
	if ~strcmpi(class(mds),'model')
		mds2=modellist;
		mds2=structtomodel(mds2,mds);
		mds=mds2;
		clear mds2;
	end
	if nargout
		varargout{1}=mds;
	else
		assignin('caller',fieldname,mds);
	end
catch me
	disp(getReport(me))
	error(['could not load model ' path]);
end
