function [Z,data_interp]=ProfileValues(md,data,xprof,yprof,resolution)
%PROFILEVALUES - compute the value of a field on a vertical profile
%
%   This routine gets the value of a given field of the model on 
%   a point given by its coordinates
%
%   Usage:
%      [z,data]=ProfileValues(md,data,xcoord,ycoord,resolution)

%Get bed and surface for each 2d point, offset to make sure that it is inside the glacier system
offset=10^-3;
bed=InterpFromMeshToMesh2d(md.mesh.elements2d,md.mesh.x2d,md.mesh.y2d,project2d(md,md.geometry.base,1),xprof,yprof)+offset;
surface=InterpFromMeshToMesh2d(md.mesh.elements2d,md.mesh.x2d,md.mesh.y2d,project2d(md,md.geometry.surface,1),xprof,yprof)-offset;

%Some useful parameters
layers=ceil(mean(md.geometry.thickness)/resolution);
Z=(bed:resolution:surface)';
X=xprof*ones(size(Z));
Y=yprof*ones(size(Z));
data_interp=InterpFromMeshToMesh3d(md.mesh.elements,md.mesh.x,md.mesh.y,md.mesh.z,data,X,Y,Z,NaN);
