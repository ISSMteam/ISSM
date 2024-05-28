function areas=GetAreasSphericalTria(index,x,y,rad_e,varargin)
%GETAREASSPHERICALTRIA - compute areas of spherical triangles 
%
%   compute areas of spherical trianguls 
%
%   Usage:
%      areas  =GetAreasSphericalTria(index,x,y,r);
%
%   Examples:
%      areas  =GetAreasSphericalTria(md.mesh.elements,md.mesh.lat,md.mesh.long,earth_radius);
%		 

%get number of elements and number of nodes
nels=size(index,1); 
nods=length(x);  

%some checks
if nargout~=1 | (nargin~=3 & nargin~=4),
	help GetAreasSphericalTria
	error('GetAreasSphericalTria error message: bad usage')
end
if (length(y)~=nods),
	error('GetAreasSphericalTria error message: x and y do not have the same length')
end
if max(index(:))>nods,
	error(['GetAreasSphericalTria error message: index should not have values above ' num2str(nods) ])
end

%initialization
areas=zeros(nels,1);
x1=x(index(:,1)); x2=x(index(:,2)); x3=x(index(:,3));
y1=y(index(:,1)); y2=y(index(:,2)); y3=y(index(:,3));

%compute the volume of each element
if nargin==4,
   % arc lengths 
	arc_12=distance(x1,y1,x2,y2).*pi./180;
	arc_23=distance(x2,y2,x3,y3).*pi./180;
	arc_31=distance(x3,y3,x1,y1).*pi./180;
	% semi perimeter 
	semi_peri=(arc_12+arc_23+arc_31)./2; 
	% spherical excess 
	excess=4*atan(sqrt(tan(semi_peri./2).*tan((semi_peri-arc_12)./2)...
		.*tan((semi_peri-arc_23)./2).*tan((semi_peri-arc_31)./2))); 
	% spherical triangle areas 
	areas=excess.*rad_e.^2; 	
end 

