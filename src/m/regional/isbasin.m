function isbasin(name)
%ISBASIN: figure out if a basin name exists.
%
%
%        Usage:  index=isbasin('jks');
%
%

%First, load basin names:
load([jplsvn '/ModelData/Names/Names.mat']);

%go through names: 
for i=1:length(names),
	if ~isempty(strfind(names{i,1},name)),
		disp(['''' names{i,1} ''' Long:' num2str(names{i,2}) ' Lat:' num2str(names{i,3}) ]);
	end
end
