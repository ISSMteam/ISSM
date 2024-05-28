function md=setmask(md,floatingicename,groundedicename,varargin)
%SETMASK - establish boundaries between grounded and floating ice.
%
%   By default, ice is considered grounded. The contour floatingicename defines nodes 
%   for which ice is floating. The contour groundedicename defines nodes inside an floatingice, 
%   that are grounded (ie: ice rises, islands, etc ...)
%   All input files are in the Argus format (extension .exp).
%
%   Usage:
%      md=setmask(md,floatingicename,groundedicename)
%
%   Examples:
%      md=setmask(md,'all','');
%      md=setmask(md,'Iceshelves.exp','Islands.exp');

%some checks on list of arguments
if ((mod(nargin,2)==0) | (nargout~=1)),
	help mask
	error('mask error message');
end

if(nargin>3)
	if(varargin(1)=='icedomain'),
		icedomainname=varargin(2);	
	else
		error('mask error message: wrong field specified. Only icedomain allowed for now.');
	end
	if ~exist(icedomainname),
        error(['setmask error message: file ' icedomainname ' not found!']);
	end
end

%Get assigned fields
x=md.mesh.x;
y=md.mesh.y;
elements=md.mesh.elements;

%Assign elementonfloatingice, elementongroundedice, vertexongroundedice and vertexonfloatingice. Only change at your own peril! This is synchronized heavily with the GroundingLineMigration module. {{{
elementonfloatingice=FlagElements(md,floatingicename);
elementongroundedice=FlagElements(md,groundedicename);

%Because groundedice nodes and elements can be included into an floatingice, we need to update. Remember, all the previous 
%arrays come from domain outlines that can intersect one another: 
elementonfloatingice=double((elementonfloatingice & ~elementongroundedice));
elementongroundedice=double(~elementonfloatingice);

%the order here is important. we choose vertexongroundedice as default on the grounding line.
vertexonfloatingice=zeros(md.mesh.numberofvertices,1);
vertexongroundedice=zeros(md.mesh.numberofvertices,1);
vertexongroundedice(md.mesh.elements(find(elementongroundedice),:))=1;
vertexonfloatingice(find(~vertexongroundedice))=1;
%}}}

%level sets
md.mask.ocean_levelset=vertexongroundedice;
md.mask.ocean_levelset(find(vertexongroundedice==0.))=-1.;

if(nargin>3)
	if(varargin(1)=='icedomain')
		md.mask.ice_levelset = 1.*ones(md.mesh.numberofvertices,1);
		%use contourtomesh to set ice values inside ice domain
		[vertexinsideicedomain,elementinsideicedomain]=ContourToMesh(elements,x,y,icedomainname,'node',1);
		pos=find(vertexinsideicedomain==1.);
		md.mask.ice_levelset(pos) = -1.;
	end
else
	md.mask.ice_levelset = -1.*ones(md.mesh.numberofvertices,1);
end


