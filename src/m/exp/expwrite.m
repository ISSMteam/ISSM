function expwrite(a,filename)
%EXPWRITE - write an Argus file from a structure given in input
%
%   This routine write an Argus file form a structure containing the fields:
%   x and y of the coordinates of the points.
%   The first argument is the structure containing the points coordinates 
%   and the second one the file to be write.
%
%   Usage:
%      expwrite(a,filename)
% 
%   Example:
%      expwrite(coordstruct,'domainoutline.exp')
%
%   See also EXPDOC, EXPREAD, EXPWRITEASVERTICES

%check input variable
if ~isstruct(a),
	error('first argument is not a structure');
end

%Add density if it's not there
if ~isfield(a,'density'),
	for n=1:length(a),
		a(n).density=1;
	end
end

fid=fopen(filename,'w');
if fid==-1,
	choice=input(['WARNING: file ' filename ' could not be created, would you like to save your exp as ./temp_expwrite.exp? (y/n)'],'s');
	if ~strcmpi(choice,'y'),
		disp('no file written... exiting');
		return
	end
	fid=fopen('./temp_expwrite.exp','w');
end
for n=1:length(a),
	if(length(a(n).x)~=length(a(n).y)),
		error('contours x and y coordinates must be of identical size');
	end

	if isfield(a,'name'),
		fprintf(fid,'%s%s\n','## Name:',a(n).name);
	else
		fprintf(fid,'%s%s\n','## Name:',filename);
	end

	fprintf(fid,'%s\n','## Icon:0');
	fprintf(fid,'%s\n','# Points Count Value');
	if isfield(a,'density'),
		if ~isempty(a(n).density),
			fprintf(fid,'%i %f\n',[length(a(n).x) a(n).density]);
		else
			fprintf(fid,'%i %f\n',[length(a(n).x) 1.]);
		end
	else
		fprintf(fid,'%i %f\n',[length(a(n).x) 1.]);
	end
	fprintf(fid,'%s\n','# X pos Y pos');
	fprintf(fid,'%10.10f %10.10f\n',[a(n).x(:) a(n).y(:)]');
	fprintf(fid,'\n');

end
fclose(fid);
