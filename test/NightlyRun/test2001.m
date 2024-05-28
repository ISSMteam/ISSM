%Test Name: SquareSheetConstrainedGia2d
%GIA test, based off of test101. Running default GIA Ivins class.
md=triangle(model(),'../Exp/Square.exp',100000.);
md=setmask(md,'','');
md=parameterize(md,'../Par/SquareSheetConstrained.par');

%GIA Ivins, 2 layer model.
md.solidearth.settings.grdmodel=2;
md.solidearth.settings.isgrd=1;

md.materials=materials('litho','ice');
md.materials.numlayers=2;
md.materials.radius =  [10 6271e3 6371e3];
md.materials.density=  [3.34e3 3.32e3];
md.materials.lame_mu=  [1.45e11         6.7e10];
md.materials.viscosity=[1e21            0];
md.initialization.sealevel=zeros(md.mesh.numberofvertices,1);
md.solidearth.settings.cross_section_shape=1;    % for square-edged x-section 
md.solidearth.settings.grdocean=0;  %do not compute sea level, only deformation
md.solidearth.settings.sealevelloading=0;  %do not compute sea level, only deformation
md.solidearth.requested_outputs={'Sealevel','BedGRD'};

%Loading history 
md.timestepping.start_time=-2400000; %4,800 kyr :: EVALUATION TIME
md.timestepping.time_step= 2400000; %2,400 kyr :: EVALUATION TIME
% to get rid of default final_time: make sure final_time>start_time
md.timestepping.final_time=2400000; %2,500 kyr
md.masstransport.spcthickness=[...
	[md.geometry.thickness; 0],...
	[md.geometry.thickness; 2400000]...
	];

%geometry at 0 initially: 
md.geometry.thickness=zeros(md.mesh.numberofvertices,1);
md.geometry.surface=zeros(md.mesh.numberofvertices,1);
md.geometry.base=zeros(md.mesh.numberofvertices,1);

%Physics: 
md.transient.issmb=0; 
md.transient.isstressbalance=0;
md.transient.isthermal=0;
md.transient.ismasstransport=1;
md.transient.isslc=1;

% Solve for GIA deflection 
md.cluster=generic('name',oshostname(),'np',3);
md.verbose.solver=0;
md=solve(md,'Transient');

%Fields and tolerances to track changes
field_names     ={'UGrd'};
field_tolerances={1e-13};
field_values={md.results.TransientSolution(2).BedGRD};
