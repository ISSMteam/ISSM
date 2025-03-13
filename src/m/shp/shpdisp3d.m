function shpdisp3d(domainoutline,varargin)
%SHPDISP - plot the contours of a domain outline file on a globe in 3d
%
%   This routine reads in a domain outline file (Shape format) and plots all the contours on a 3D rendition of the earth.
%
%   Usage:
%      shpdisp3d(domainoutline,varargin)
%      shpdisp3d(filenamei,'figure',1,'style',stylei,'linewidth',linewidthi);
%
%   Example:
%      shpdisp3d('Domain.shp','figure',1,'style','--r','linewidthi',2);
%
%   See also SHPREAD, SHPDOC, SHPDISP

%recover options
options=pairoptions(varargin{:});

%parse input:
figurenumber=getfieldvalue(options,'figure',1);
color=getfieldvalue(options,'color','r');
linewidth=getfieldvalue(options,'linewidth',1);
unitmultiplier=getfieldvalue(options,'unit',1);
epsg=getfieldvalue(options,'epsg',4326);
radius=getfieldvalue(options,'radius',6371012);
aboveground=getfieldvalue(options,'aboveground',1000)

%read domain:
domain=shpread(domainoutline);

if epsg~=4326,
	%transform to lat,long:
	for i=1:length(domain),
		[x,y] = CoordTransform(domain(i).x,domain(i).y,'EPSG:4326',sprintf('EPSG:%i',epsg));
		domain(i).x=x; domain(i).y=y;
	end
end

for i=1:length(domain),

	%make sure lat,long are what they are supposed to be: 
	if any(domain(i).x>90 | domain(i).x<-90), 
		long=domain(i).x; lat=domain(i).y;
	else
		long=domain(i).y; lat=domain(i).x;
	end

	%project on x,y,z reference frame.
	[x,y,z]=AboveGround(lat,long,radius,aboveground);
	hold on, p=plot3(x,y,z,'k-'); 
	set(p,'Color',color);
	set(p,'LineWidth',linewidth);
end
