%Test Name: ElasticLoveNumbers
%Forward Love number solution for a viscoelastic earth, model M3-L70-V01 from
%Spada, G., Barletta, V. R., Klemann, V., Riva, R. E. M., Martinec, Z.,
%Gasperini, P., Lund, B., Wolf, D., Vermeersen, L. L. A. and King, M. A.
%(2011), A benchmark study for glacial isostatic adjustment codes. Geophysical
%Journal International, 185: 106--132. doi:10.1111/j.1365-246X.2011.04952.x

md=model();
md.cluster=generic('name',oshostname(),'np',3);

md.materials=materials('litho');
md.miscellaneous.name='test2070';
md.groundingline.migration='None';

md.verbose=verbose('all');
md.verbose=verbose('1111111111111111');
cst=365.25*24*3600*1000;

md.materials.numlayers=6;
md.materials.radius =  [10 1222.5 3.4800e+03   5.7010e+03   5.9510e+03   6.3010e+03   6.3710e+03]'*1e3;
md.materials.density=  [1.0750e4 1.0750e+04   4.9780e+03   3.8710e+03   3.4380e+03   3.0370e+03]';
md.materials.lame_mu=  [1e-5         0   2.2834e+00   1.0549e+00   7.0363e-01   5.0605e-01]'*1e11;
md.materials.viscosity=[0            0   2.0000e+00   1.0000e+00   1.0000e+00   1.0000e+25]'*1e21;
md.materials.lame_lambda=md.materials.lame_mu*0+5e17;
md.materials.issolid=[1 0 1 1 1 1]';
md.materials.rheologymodel=zeros(md.materials.numlayers,1);

md.love.allow_layer_deletion=1;
md.love.frequencies=[0];
md.love.nfreq=length(md.love.frequencies);
md.love.sh_nmin=1;
md.love.sh_nmax=10000;
md.love.underflow_tol=1e-20;
md.love.pw_threshold=1e-3;
md.love.Gravitational_Constant=6.6732e-11;
md.love.min_integration_steps=100;
md.love.istemporal=0;
%md.love.time=(linspace(1/12,10, 10*12))'*cst/1e3;

md=solve(md,'lv');

h=md.results.LoveSolution.LoveHf;
l=md.results.LoveSolution.LoveLf;
k=md.results.LoveSolution.LoveKf;

md.love.forcing_type=9;
md.love.sh_nmin=2;
md=solve(md,'lv');

th=md.results.LoveSolution.LoveHf;
tl=md.results.LoveSolution.LoveLf;
tk=md.results.LoveSolution.LoveKf;

%Fields and tolerances to track changes

field_names     ={'LoveH_loading_elastic','LoveK_loading_elastic','LoveL_loading_elastic','LoveH_tidal_elastic','LoveK_tidal_elastic','LoveL_tidal_elastic'};
field_tolerances={2.0e-8,2.0e-8,2.0e-8,2.0e-8,2.0e-8,2.0e-8};
field_values={h,k,l,th,tk,tl};

