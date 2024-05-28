function [intData, meanData, areas] = integrateOverDomain(md, data, masked, weights)
% integrateOverDomain - integrating data over the whole domain
%
%   intData: integral of the data over each element
%   meanData: intData/areas
%   areas: areas of the domain
if nargin < 4
	weights = ones(size(data));
	if nargin<3
		masked = logical(zeros(size(data)));
	end
end

masked = masked | isnan(data) | isnan(weights);
% Set the area with masked=1 to nan
data(masked) = nan;
weights(masked) =nan;


% get the mesh
elements=md.mesh.elements;
x=md.mesh.x;
y=md.mesh.y;

%compute areas;
eleAreas=GetAreas(elements,x,y);

% integrate nodal data to element
eleData = 1/3*eleAreas.*(data(elements(:,1),:).*weights(elements(:,1),:) + data(elements(:,2),:).*weights(elements(:,2),:) + data(elements(:,3),:).*weights(elements(:,3),:));
eleAreas = 1/3*eleAreas.*(weights(elements(:,1),:)+weights(elements(:,2),:)+weights(elements(:,3),:));

intData = sum(eleData, 1, 'omitnan');
areas = sum(eleAreas, 1, 'omitnan');
meanData = intData ./ areas;
