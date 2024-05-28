%Test Name: GiaIvinsBenchmarksAB2dA
% Benchmark experiments (Figure A2a Ivins and James, 1999, Geophys. J. Int.) 
md=triangle(model(),'../Exp/RoundFrontEISMINT.exp',200000.);
md=setmask(md,'','');
md=parameterize(md,'../Par/GiaIvinsBenchmarksAB.par');

%GIA Ivins, 2 layer model.
md.solidearth.settings.grdmodel=2;
md.solidearth.settings.isgrd=1;
md.initialization.sealevel=zeros(md.mesh.numberofvertices,1);

%% indicate what you want to compute 
md.solidearth.settings.cross_section_shape=1;    % for square-edged x-section 
md.solidearth.settings.grdocean=0;  %do not compute sea level, only deformation
md.solidearth.settings.sealevelloading=0;  %do not compute sea level, only deformation

% evaluation time (termed start_time) 

md.timestepping.time_step=2002100; % after 2 kyr of deglaciation 
md.timestepping.start_time=-md.timestepping.time_step; % need one time step before t=0 to generate a thickness change at t=0
md.timestepping.final_time=2500000; % 2,500 kyr

%% define loading history {{{ 
md.masstransport.spcthickness=[...
	[md.geometry.thickness*0.0; 0.0],...
	[md.geometry.thickness; 1000],...
	[md.geometry.thickness; 2000000],...
	[md.geometry.thickness*0.0; 2000100],...
	[md.geometry.thickness*0.0; md.timestepping.start_time+2*md.timestepping.time_step],...
	];
% }}} 

md.geometry.bed=zeros(md.mesh.numberofvertices,1);

% find out elements that have zero loads throughout the loading history.
pos = find(sum(abs(md.geometry.thickness(1:end-1,:)),2)==0);
md.mask.ice_levelset(pos)=1; % no ice

%Physics: 
md.transient.issmb=0; 
md.transient.isstressbalance=0;
md.transient.isthermal=0;
md.transient.ismasstransport=1;
md.transient.isslc=1;

md.cluster=generic('name',oshostname(),'np',3);
md.verbose=verbose('1111111');
md.solidearth.requested_outputs={'Sealevel','BedGRD'};

%% solve for GIA deflection 
md=solve(md,'Transient');

%Test Name: GiaIvinsBenchmarksAB2dA1
U_AB2dA1 = md.results.TransientSolution.BedGRD; 
%URate_AB2dA1 = md.results.TransientSolution.UGiaRate; 

%Test Name: GiaIvinsBenchmarksAB2dA2
%% different evaluation time. {{{ 
md.timestepping.time_step=2005100;% after 5 kyr of deglaciation
md.timestepping.start_time=-md.timestepping.time_step; % need one time step before t=0 to generate a thickness change at t=0
md.masstransport.spcthickness(end,end) = md.timestepping.start_time+2*md.timestepping.time_step;

md=solve(md,'Transient');

U_AB2dA2 = md.results.TransientSolution.BedGRD; 
%URate_AB2dA2 = md.results.TransientSolution.BedGRDRate; 
% }}} 

%Test Name: GiaIvinsBenchmarksAB2dA3
%% different evaluation time. {{{ 
md.timestepping.time_step=2010100;% after 10 kyr of deglaciation
md.timestepping.start_time=-md.timestepping.time_step; % need one time step before t=0 to generate a thickness change at t=0
md.masstransport.spcthickness(end,end) = md.timestepping.start_time+2*md.timestepping.time_step; 

md=solve(md,'Transient');

U_AB2dA3 = md.results.TransientSolution.BedGRD; 
%URate_AB2dA3 = md.results.TransientSolution.UGiaRate; 
% }}} 

%Fields and tolerances to track changes
field_names     ={'U_AB2dA1','U_AB2dA2','U_AB2dA3'};
field_tolerances={1e-13,1e-13,1e-13};
field_values={U_AB2dA1,U_AB2dA2,U_AB2dA3};
% field_names     ={'U_AB2dA1','URate_AB2dA1','U_AB2dA2','URate_AB2dA2','U_AB2dA3','URate_AB2dA3'};
% field_tolerances={1e-13,1e-13,1e-13,1e-13,1e-13,1e-13};
% field_values={U_AB2dA1,URate_AB2dA1,U_AB2dA2,URate_AB2dA2,U_AB2dA3,URate_AB2dA3};

