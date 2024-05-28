%Test Name: SquareSheetShelfGroundingLine2dAggressive. From test424, with sea level increasing.
md=triangle(model(),'../Exp/Square.exp',150000.);
md=setmask(md,'../Exp/SquareShelf.exp','');
md=parameterize(md,'../Par/SquareSheetShelf.par');
md=setflowequation(md,'SSA','all');
md.initialization.vx(:)=0.;
md.initialization.vy(:)=0.;
md.smb.mass_balance(:)=0.;

md.geometry.base=-700.-abs(md.mesh.y-500000.)/1000.;
md.geometry.bed =-700.-abs(md.mesh.y-500000.)/1000.;
md.geometry.thickness(:)=1000.;
md.geometry.surface=md.geometry.base+md.geometry.thickness;

md.transient.isstressbalance=0;
md.transient.isgroundingline=1;
md.transient.isthermal=0;
md.groundingline.migration='AggressiveMigration';
md.transient.requested_outputs={'IceVolume','IceVolumeAboveFloatation','Sealevel'};

md.timestepping.time_step=.1;
md.solidearth.initialsealevel=newforcing(md.timestepping.start_time, md.timestepping.final_time, md.timestepping.time_step,-200,200,md.mesh.numberofvertices);

md.cluster=generic('name',oshostname(),'np',3);
md=solve(md,'Transient');

%we are checking that the grounding line position is near the theorical one, which is the 0 contour level 
%of surface - sealevel - (1-di)* thickness; 

nsteps=length(md.results.TransientSolution);
field_names={};
field_tolerances={};
field_values={};
for i=1:nsteps,
	field_names{end+1}=sprintf('Time-%g-yr-ice_levelset-S-sl-(1-di)*H',md.results.TransientSolution(i).time);
	field_tolerances{end+1}=1e-12;
	field_values{end+1}= md.results.TransientSolution(i).MaskOceanLevelset- (md.geometry.surface-md.results.TransientSolution(i).Sealevel-(1-md.materials.rho_ice/md.materials.rho_water)*md.geometry.thickness);
end
