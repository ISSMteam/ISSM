function Struct=meshread(filename)

%some checks
if ~exist(filename),
	error(['meshread error message: file ' filename ' not found!']);
end

fid=fopen(filename,'r');

while (~feof(fid)),

	A=fscanf(fid,'%s',1);

	if strcmp(A,'MeshVersionFormatted');
		Struct.Version=fscanf(fid,'%s',1);

	elseif strcmp(A,'Dimension'),
		Struct.Dimension=fscanf(fid,'%i',1);

	elseif strcmp(A,'Vertices'),
		Struct.nods=fscanf(fid,'%i',1);
		A=fscanf(fid,'%f %f %f',[3 Struct.nods]);
		Struct.x=A(1,:)';
		Struct.y=A(2,:)';

	elseif strcmp(A,'Triangles'),
		Struct.nels=fscanf(fid,'%i',1);
		A=fscanf(fid,'%i %i %i',[4 Struct.nels]);
		Struct.index=A(1:3,:)';

	elseif strcmp(A,'Quadrilaterals'),
		Struct.nels=fscanf(fid,'%i',1);
		A=fscanf(fid,'%i %i %i %i',[5 Struct.nels]);
		Struct.index=A(1:4,:)';
	else
		%do nothing

	end
end

fclose(fid);
