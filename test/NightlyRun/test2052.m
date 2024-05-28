%Test Name: GiaIvinsBenchmarksAB2dC
% Benchmark experiments (Figure A2a Ivins and James, 1999, Geophys. J. Int.) 
md=triangle(model(),'../Exp/RoundFrontEISMINT.exp',200000.);
md=setmask(md,'','');
md=parameterize(md,'../Par/GiaIvinsBenchmarksCD.par');

%% indicate what you want to compute 
md.gia.cross_section_shape=1;    % for square-edged x-section 

% evaluation time (termed start_time) 
md.timestepping.start_time=0.3;     % for t \approx 0 kyr: to get elastic response!
md.timestepping.final_time=2500000; % 2,500 kyr

%% define loading history {{{ 
md.geometry.thickness=[...
	[md.geometry.thickness*0.0; 0.0],...
	[md.geometry.thickness/2.0; 0.1],...
	[md.geometry.thickness; 0.2],...
	[md.geometry.thickness; md.timestepping.start_time],...
	];
% }}} 

% find out elements that have zero loads throughout the loading history.
pos = find(sum(abs(md.geometry.thickness(1:end-1,:)),2)==0);
md.mask.ice_levelset(pos)=1; % no ice

md.cluster=generic('name',oshostname(),'np',3);
md.verbose=verbose('1111111');

%% solve for GIA deflection 
md=solve(md,'Gia');

%Test Name: GiaIvinsBenchmarksAB2dC1
U_AB2dC1 = md.results.GiaSolution.UGia; 
URate_AB2dC1 = md.results.GiaSolution.UGiaRate;

%Test Name: GiaIvinsBenchmarksAB2dC2
%% different evaluation time. {{{ 
md.timestepping.start_time=1000.3;  % for t \approx 1 kyr 
md.geometry.thickness(end,end) = md.timestepping.start_time; 

md=solve(md,'Gia');

U_AB2dC2 = md.results.GiaSolution.UGia; 
URate_AB2dC2 = md.results.GiaSolution.UGiaRate; 
% }}} 

%Test Name: GiaIvinsBenchmarksAB2dC3
%% different evaluation time. {{{ 
md.timestepping.start_time=2400000; % for t \approx \infty 
md.geometry.thickness(end,end) = md.timestepping.start_time; 

md=solve(md,'Gia');

U_AB2dC3 = md.results.GiaSolution.UGia; 
URate_AB2dC3 = md.results.GiaSolution.UGiaRate; 
% }}} 

%Fields and tolerances to track changes
field_names     ={'U_AB2dC1','URate_AB2dC1','U_AB2dC2','URate_AB2dC2','U_AB2dC3','URate_AB2dC3'};
field_tolerances={1e-13,1e-13,1e-13,1e-13,1e-13,1e-13};
field_values={U_AB2dC1,URate_AB2dC1,U_AB2dC2,URate_AB2dC2,U_AB2dC3,URate_AB2dC3}; 

