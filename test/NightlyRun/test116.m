%Test Name: SquareShelfConstrainedBalThic2d
md=triangle(model(),'../Exp/Square.exp',150000.);
md=setmask(md,'all','');
md=parameterize(md,'../Par/SquareShelfConstrained.par');
%Add boundary conditions on thickness on the border
pos=find(md.mesh.vertexonboundary);
md.balancethickness.spcthickness(pos)=md.geometry.thickness(pos);
md=setflowequation(md,'SSA','all');
md.cluster=generic('name',oshostname(),'np',3);
md=solve(md,'Balancethickness');

%Fields and tolerances to track changes
field_names     ={'Thickness'};
field_tolerances={1e-13};
field_values={...
	(md.results.BalancethicknessSolution.Thickness),...
	};
