%Test Name: SquareSheetShelfStressFSEstar
md=triangle(model(),'../Exp/Square.exp',180000.);
md=setmask(md,'../Exp/SquareShelf.exp','');
md=parameterize(md,'../Par/SquareSheetShelf.par');
md=extrude(md,3,1.);
md.materials = matestar(md.materials);
md.materials.rheology_B = 3.15e8*ones(md.mesh.numberofvertices,1);
md.materials.rheology_Ec=ones(md.mesh.numberofvertices,1);
md.materials.rheology_Es=3*ones(md.mesh.numberofvertices,1);
md.cluster=generic('name',oshostname(),'np',3);

%Go solve
field_names={};
field_tolerances={};
field_values={};
%md.initialization.pressure=md.constants.g*md.materials.rho_ice*(md.geometry.surface-md.mesh.y);
for i={'SSA','HO','FS'},
	disp(' ');
	disp(['====== Testing Estar with ' i{1} ' =====']);
	md=setflowequation(md,i{1},'all');
	md=solve(md,'Stressbalance');
	field_names     ={field_names{:},['Vx' i{1}],['Vy' i{1}],['Vz' i{1}],['Vel' i{1}],['Pressure' i{1}]};
	field_tolerances={field_tolerances{:},7e-06,2e-05,2e-06,5e-06,8e-07};
	field_values={field_values{:},...
		(md.results.StressbalanceSolution.Vx),...
		(md.results.StressbalanceSolution.Vy),...
		(md.results.StressbalanceSolution.Vz),...
		(md.results.StressbalanceSolution.Vel),...
		(md.results.StressbalanceSolution.Pressure),...
		};
end
