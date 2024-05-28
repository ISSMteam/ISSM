function [alpha,beta,theta,ex,sigxx]=thomasparams(md,varargin)
%THOMASPARAMS - compute Thomas' geometric parameters for an ice shelf 
%
%	 This routine computes geometric parameters representing ratios between
%	 components of the horizontal strain rate tensor for an ice shelf, as
%	 originally developed in Thomas (1973).  The model must contain computed
%	 strain rates, either from observed or modeled ice velocities.
%
%   Available options:
%	 -'eq'			: analytical equation to use in the calculation.  Must be one of:
%				'Thomas' for a 2D ice shelf, taking into account full strain rate
%					tensor (default)
%				'Weertman1D' for a confined ice shelf free to flow in one direction
%				'Weertman2D' for an unconfined ice shelf free to spread in any direction
%
%	 -'smoothing'	: an integer smoothing parameter for the averaging function
%						(default 0) Type 'help averaging' for more information on its usage.
%
%	 -'coordsys'	: coordinate system for calculating the strain rate
%						components. Must be one of:
%				'longitudinal': x axis aligned along a flowline at every point (default)
%				'principal': x axis aligned along maximum principal strain rate
%					at every point
%				'xy': x and y axes same as in polar stereographic projection 
%
%   Return values: 
%
%		'alpha' which is the ratio e_yy/e_xx between components of the strain
%		rate tensor
%
%		'beta' which is the ratio e_xy/e_xx between components of the strain rate
%		tensor
%
%		'theta' which is a combination of alpha and beta arising from the form of
%		the equivalent stress
%
%		'exx' is the strain rate along a coordinate system defined by 'coordsys' 
%
%		'sigxx' is the deviatoric stress along a coordinate system defined by 'coordsys' 
%
%   Usage: [alpha,beta,theta,exx,sigxx]=ThomasParams(md,options)
%
%   Example: [alpha,beta,theta,exx,sigxx]=ThomasParams(md,'eq','Thomas','smoothing',2,'coordsys','longitudinal')

%some checks
if (nargin<4)
	help ThomasParams
	error('bad usage');
end
if isempty(fieldnames(md.results)),
	error(['md.results.strainrate is not present.  Calculate using md=mechanicalproperties(md,vx,vy)'])
end
if dimension(md.mesh)~=2,
	error('only 2d model supported currently');
end
if any(md.flowequation.element_equation~=2),
	disp('Warning: the model has some non SSA elements. These will be treated like SSA elements');
end

% process options
options = pairoptions(varargin{:});
eq = getfieldvalue(options,'eq','Thomas');
smoothing = getfieldvalue(options,'smoothing',0);
coordsys = getfieldvalue(options,'coordsys','longitudinal');

% average element strain rates onto vertices
e1=averaging(md,md.results.strainrate.principalvalue1,smoothing)/md.constants.yts; % convert to s^-1
e2=averaging(md,md.results.strainrate.principalvalue2,smoothing)/md.constants.yts;
exx=averaging(md,md.results.strainrate.xx,smoothing)/md.constants.yts;
eyy=averaging(md,md.results.strainrate.yy,smoothing)/md.constants.yts;
exy=averaging(md,md.results.strainrate.xy,smoothing)/md.constants.yts;

% checks: any of e1 or e2 equal to zero?
pos=find(e1==0);
if any(pos==1)
	disp('WARNING: first principal strain rate equal to zero.  Value set to 1e-13 s^-1');
	e1(pos)=1e-13;
end
pos=find(e2==0);
if any(pos==1)
	disp('WARNING: second principal strain rate equal to zero.  Value set to 1e-13 s^-1');
	e2(pos)=1e-13;
end

% rheology
n=averaging(md,md.materials.rheology_n,0);
B=md.materials.rheology_B;

switch coordsys
	case 'principal'
		b=zeros(md.mesh.numberofvertices,1);
		ex=e1;
		a=e2./e1;
		pos=find(e1<0 & e2>0); % longitudinal compression and lateral tension
		a(pos)=e1(pos)./e2(pos);
		ex(pos)=e2(pos);
		pos2=find(e1<0 & e2<0 & abs(e1)<abs(e2)); % lateral and longitudinal compression
		a(pos2)=e1(pos2)./e2(pos2);
		ex(pos2)=e2(pos2);
		pos3=find(e1>0 & e2>0 & abs(e1)<abs(e2)); % lateral and longitudinal tension
		a(pos3)=e1(pos3)./e2(pos3);
		ex(pos3)=e2(pos3);
		id=find(e1<0 & e2<0);
		a(id)=-a(id); % where both strain rates are compressive, enforce negative alpha
		sigxx=(abs(ex)./((1+a+a.^2).^((n-1)/2))).^(1./n).*B;
		%mask=ismember(1:md.mesh.numberofvertices,id);
		%plotmodel(md,'data',ex,'mask',mask)
	case 'xy'
		ex=exx;
		a=eyy./exx;
		b=exy./exx;
	case 'longitudinal'
		% using longitudinal strain rates defined by observed velocity vector
		velangle=atan(md.initialization.vy./md.initialization.vx);
		pos=find(md.initialization.vx==0);
		velangle(pos)=pi/2;
		ex=0.5*(exx+eyy)+0.5*(exx-eyy).*cos(2*velangle)+exy.*sin(2*velangle);
		ey=exx+eyy-ex; % trace of strain rate tensor is invariant
		exy=-0.5*(exx-eyy).*sin(2*velangle)+exy.*cos(2*velangle);
		a=ey./ex;
		b=exy./ex;
		%pos=find(ex<0 & ey<0);
		%a(pos)=-a(pos);
		%sigxx=(abs(ex)./((1+a+a.^2+b.^2).^((n-1)/2))).^(1./n).*B;
		sigxx=abs(ex).^(1./n-1).*ex./((1+a+a.^2+b.^2).^((n-1)./(2*n))).*B;
	otherwise
		error('argument passed to "coordsys" not valid');
end

% a < -1 in areas of strong lateral compression or longitudinal compression and
% theta flips sign at a = -2
pos=find(abs((abs(a)-2))<1e-3);
if length(pos)>0,
	disp(['Warning: ', num2str(length(pos)), ' vertices have alpha within 1e-3 of -2'])
end
a(pos)=-2+1e-3;

switch eq
	case 'Weertman1D'
		theta=1./8;
		a=zeros(md.mesh.numberofvertices,1);
	case 'Weertman2D'
		theta=1./9;
		a=ones(md.mesh.numberofvertices,1);
	case 'Thomas'
		theta=((1+a+a.^2+b.^2).^((n-1)/2))./(abs(2+a).^n);
	otherwise
		error('argument passed to "eq" not valid.  Type "help zinv" for usage');
end
alpha=a;
beta=b;
