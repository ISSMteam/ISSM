function [data_interp]=PointValues(md,data,xpoint,ypoint)
%POINTVALUES - compute the value of a field on a single point
%
%   This routine gets the value of a given field of the model on points
%
%   Usage:
%      [z,data]=PointValues(md,data,X,Y,Z)

%Get bed and surface for each 2d point, offset to make sure that it is inside the glacier system
%offset=10^-3;
%bed=InterpFromMeshToMesh2d(md.mesh.elements2d,md.mesh.x2d,md.mesh.y2d,project2d(md,md.geometry.bed,1),xpoint,ypoint)+offset;
%surface=InterpFromMeshToMesh2d(md.mesh.elements2d,md.mesh.x2d,md.mesh.y2d,project2d(md,md.geometry.surface,1),xpoint,ypoint)-offset;

data_interp=InterpFromMeshToMesh2d(md.mesh.elements,md.mesh.x,md.mesh.y,data,xpoint,ypoint);
