%ga: grounded area
%iv: ice volume
%ivaf: ice volume above floatation
%GL_y : grounding line position @ y (y comes in m)
%nelem : number of elements
% usage:
%
% Default: y=40km, i0=1
% [ga iv ivaf GL_y IF_y nelem t] = ice_evolution(md);
%
% Default: y=40km
% [ga iv ivaf GL_y IF_y nelem t] = ice_evolution(md,i0);
%
% Use this for y the borders (y=0 or y=ymax)
% [ga iv ivaf GL_y IF_y nelem t] = ice_evolution(md,i0,y);
%
%

function [ga iv ivaf GL_y IF_y nelem t] = ice_evolution(varargin),

	ga			= [];
	iv			= [];
	ivaf		= [];
	GL_y		= [];
	IF_y		= [];
	nelem		= [];
	t			= [];

	if(nargin==0)
		error('it is necessary the model!')
	elseif(nargin==1)
		% Defoult is y=40km
		i0	= 1;
		y0	= 35000;
		y1 = 45000;
		dy = 100;
		y  = 40000;
	elseif(nargin==2)
		i0=varargin{2};
		% Defoult is y=40km
		y0	= 35000;
		y1 = 45000;
		dy = 100;
		y  = 40000;
	elseif(nargin==3)
		i0 = varargin{2};
		y  = varargin{3};
		dy = 10;
		y0	= y-dy;
		y1 = y+dy;
	else
		error('number of inputs is more than 3');
	end
	%set the model
	md			=varargin{1};
	nsteps	= length(md.results.TransientSolution);


	for i=i0:nsteps,
		ga(i)			= md.results.TransientSolution(i).GroundedArea;
		iv(i)			= md.results.TransientSolution(i).IceVolume;
		ivaf(i)		= md.results.TransientSolution(i).IceVolumeAboveFloatation;
		if(isfield(md.results.TransientSolution,'MeshElements'))
			nelem(i)		= size(md.results.TransientSolution(i).MeshElements,1);
		else
			nelem(i)		= md.mesh.numberofelements;
		end
		t(i)			= md.results.TransientSolution(i).time;	
		%find GL position between y0 and y1 {{{
		[glx gly]	= gl_position(md,i,0);
		pos			= find(gly<y1 & gly>y0);
		x				= gly(pos);
		v				= glx(pos);
		if(length(pos)==0)
			error('pos is null')
		elseif(length(pos)==1)
			%this should be used for y=0 or y=ymax
			GL_y(i)	= v;
		else
			%this should be used when y is inside the domain; so, use linear interpolation
			xq			= [y0:dy:y1];
			vq			= interp1(x,v,xq,'linear');
			pos		= find(xq==y);
			if(pos)
				GL_y(i)	= vq(pos);
			else
				error('pos is null')
			end
		end
		%}}}	
		%find IF position between y0 and y1 {{{
		[ifx ify]	= if_position(md,i,0);
		if(length(ifx)==0)
			continue
		end
		pos			= find(ify<y1 & ify>y0);
		x				= ify(pos);
		v				= ifx(pos);
		if(length(pos)==0)
			error('pos is null')
		elseif(length(pos)==1)
			%this should be used for y=0 or y=ymax
			IF_y(i)	= v;
		else
			%this should be used when y is inside the domain; so, use linear interpolation
			xq			= [y0:dy:y1];
			vq			= interp1(x,v,xq,'linear');
			pos		= find(xq==y);
			if(pos)
				IF_y(i)	= vq(pos);
			else
				error('pos is null')
			end
		end
		%}}}
	end

end
