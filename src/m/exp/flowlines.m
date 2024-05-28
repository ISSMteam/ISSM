function flowpath=flowlines(index,x,y,u,v,x0,y0,varargin)
%FLOWLINES - compute flowlines from a given set of seed points
%
%   Usage:
%      flowpath=flowlines(index,x,y,u,v,x0,y0,varargin)
%
%   the velocity field is given by the couple (u,v) and the coordinates
%   of the seed points are (x0,y0). One can use one or several seed 
%   points
%
%   Example:
%      flowpath=flowlines(md.mesh.elements,md.mesh.x,md.mesh.y,md.initialization.vx,md.initialization.vy,x0,y0)
%
%   Options:
%      - 'maxiter':   how many steps upstream and downstream of the seed points (default: 200)
%      - 'precision': division of each segment (higer precision increases number of segments, default: 1)
%      - 'downstream':flow line upstream of the seed points (default: 1)
%      - 'upstream':  flow line upstream of the seed points (default: 1)

%check input
if (length(x)~=length(y) | length(x)~=length(u) | length(x)~=length(v)),
	error('flowlines error message: x,y,u and v must have the same length');
end
if length(x)<3,
	error('flowlines error message: at least one element is required');
end
if length(x0)~=length(y0),
	error('flowlines error message: x0 and y0 do not have the same length');
end

%process options
options    = pairoptions(varargin{:});
maxiter    = getfieldvalue(options,'maxiter',200);
precision  = getfieldvalue(options,'precision',1);
downstream = getfieldvalue(options,'downstream',1);
upstream   = getfieldvalue(options,'upstream',1);

%Create triangulation once for all and check seed points
trep = triangulation(index,x,y);
tria = pointLocation(trep,[x0 y0]);
pos=find(isnan(tria));
x0(pos)=[];
y0(pos)=[];

%initialize other variables
N=length(x0);
X=x0; Y=y0;
flowpath=struct('x',cell(N,1),'y',cell(N,1),'name','','density',1);
for i=1:N,
	flowpath(i).x=x0(i);
	flowpath(i).y=y0(i);
end
done=zeros(N,1);

%get avegared length of each element
length_tria=1/3*(sqrt( (x(index(:,1))-x(index(:,2))).^2+(y(index(:,1))-y(index(:,2))).^2 )+...
	sqrt((x(index(:,1))-x(index(:,3))).^2+(y(index(:,1))-y(index(:,3))).^2 )+...
	sqrt((x(index(:,2))-x(index(:,3))).^2+(y(index(:,2))-y(index(:,3))).^2 ));

%take velocity for each element
u=u(index)*[1;1;1]/3;
v=v(index)*[1;1;1]/3;

if downstream,
	%initialization:
	counter=1;

	while any(~done) 

		%find current triangle
		queue=find(~done);
		tria = pointLocation(trep,[X(queue),Y(queue)]);

		%check that the point is actually inside a triangle of the mesh
		listnan=find(isnan(tria));
		for i=1:length(listnan)
			%remove the last point
			flowpath(queue(listnan(i))).x(end)=[];
			flowpath(queue(listnan(i))).y(end)=[];
			done(queue(listnan(i)))=1;
		end
		tria(listnan)=[]; 
		queue(listnan)=[];

		if isempty(tria),
			break;
		end

		%velocity of the current triangle and norm it
		ut=u(tria); vt=v(tria); normv=max(eps,sqrt(ut.^2+vt.^2));
		ut=ut./normv;vt=vt./normv;

		%check counter
		if counter>maxiter
			disp(['Maximum number of iterations (' num2str(maxiter) ') reached while going forward'])
			break
		end
		counter=counter+1;

		%remove stagnant point
		done(queue(find(ut==0 & vt==0)))=1;

		%build next point
		for i=1:length(queue)
			X(queue(i))=flowpath(queue(i)).x(end)+ut(i)*length_tria(tria(i))/precision;
			Y(queue(i))=flowpath(queue(i)).y(end)+vt(i)*length_tria(tria(i))/precision;
			flowpath(queue(i)).x=[flowpath(queue(i)).x;flowpath(queue(i)).x(end)+ut(i)*length_tria(tria(i))/precision];
			flowpath(queue(i)).y=[flowpath(queue(i)).y;flowpath(queue(i)).y(end)+vt(i)*length_tria(tria(i))/precision];
		end
	end
end

%same process but reverse (vel=-vel) to have a vcomplete flow line
if upstream,
	queue=[];
	counter=1;
	X=x0; Y=y0;
	done=zeros(N,1);

	while any(~done) 

		%find current triangle
		queue=find(~done);
		tria = pointLocation(trep,[X(queue),Y(queue)]);

		%check that the point is actually inside a triangle of the mesh
		listnan=find(isnan(tria));
		for i=1:length(listnan)
			%remove the last point
			flowpath(queue(listnan(i))).x(1)=[];
			flowpath(queue(listnan(i))).y(1)=[];
			done(queue(listnan(i)))=1;
		end
		tria(listnan)=[]; 
		queue(listnan)=[];

		if isempty(tria),
			break;
		end

		%velocity of the current triangle and norm it
		ut=-u(tria); vt=-v(tria); normv=max(eps,sqrt(ut.^2+vt.^2));
		ut=ut./normv;vt=vt./normv;

		%check counter
		if counter>maxiter
			disp(['Maximum number of iterations (' num2str(maxiter) ') reached while going backward'])
			break
		end
		counter=counter+1;

		%remove stagnant point
		done(queue(find(ut==0 & vt==0)))=1;

		%build next point
		for i=1:length(queue)
			X(queue(i))=flowpath(queue(i)).x(1)+ut(i)*length_tria(tria(i))/precision;
			Y(queue(i))=flowpath(queue(i)).y(1)+vt(i)*length_tria(tria(i))/precision;
			flowpath(queue(i)).x=[flowpath(queue(i)).x(1)+ut(i)*length_tria(tria(i))/precision; flowpath(queue(i)).x];
			flowpath(queue(i)).y=[flowpath(queue(i)).y(1)+vt(i)*length_tria(tria(i))/precision; flowpath(queue(i)).y];
		end
	end
end

%EXP compatibility (add name)
for i=1:length(queue)
	flowpath(queue(i)).name=['flowline' num2str(i)];
end
