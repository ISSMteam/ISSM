%Test Name: SquareBamgMesh

%Simple mesh
md=bamg(model(),'domain','../Exp/Square.exp','hmax',100000.);
x1=md.mesh.x;
y1=md.mesh.y;

%hVertices
md=bamg(model(),'domain','../Exp/Square.exp','hmax',300000.,'hVertices',[10000. 100000. 400000. 100000.]');
x2=md.mesh.x;
y2=md.mesh.y;

%big mesh
t0=clock;
md=bamg(model(),'domain','../Exp/Square.exp','hmax',3000.);
nbelements=md.mesh.numberofelements;
if nbelements>267895-50 & nbelements<267895+50
	nbewithinrange = 1.;
else
	nbewithinrange = 0.;
end
elapsedtime=etime(clock,t0);

%Fields and tolerances to track changes
field_names     ={'x1','y1','x2','y2','nbelements','elapsed time'};
field_tolerances={2e-9,2e-9,1e-13,1e-13,1e-13,8.5};
field_values={...
	x1, y1,...
	x2, y2,...
	nbewithinrange,elapsedtime...
	};
