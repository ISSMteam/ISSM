%Test Name: ElasticLoveNumbers_HighlyLayeredEarth
%Forward Love number solution for a viscoelastic earth, model M3-L70-V01 from
%Spada, G., Barletta, V. R., Klemann, V., Riva, R. E. M., Martinec, Z.,
%Gasperini, P., Lund, B., Wolf, D., Vermeersen, L. L. A. and King, M. A.
%(2011), A benchmark study for glacial isostatic adjustment codes. Geophysical
%Journal International, 185: 106--132. doi:10.1111/j.1365-246X.2011.04952.x

md=model();
md.cluster=generic('name',oshostname(),'np',3);

md.materials=materials('litho');
md.miscellaneous.name='test2072';
md.groundingline.migration='None';

md.verbose=verbose('all');
md.verbose=verbose('1111111111111111');
cst=365.25*24*3600*1000;


%Model VSS96 from Vermeersen, L.L.A., Sabadini, R. & Spada, G., 1996a. Analytical visco-elastic relaxation models, Geophys. Res. Lett., 23, 697–700.
md.materials.radius=[10, 1222.5, 3480., 3600., 3630.5, 3700., 3900., 4000., 4200., 4300., 4500., 4600., 4800., 4900., 5100., 5200., 5400., 5500., 5600.5, 5650., 5701., 5736., 5771.5, 5821., 5951., 5970.5, 6016., 6061., 6150.5, 6151.5, 6251., 6371.]'*1e3;
md.materials.lame_mu=[1e-5, 0., 2.933, 2.8990002, 2.8550003, 2.7340002, 2.675, 2.559, 2.502, 2.388, 2.331, 2.215, 2.157, 2.039, 1.979, 1.8560001, 1.794, 1.73, 1.639, 1.2390001, 1.224, 1.21, 1.128, 0.97700006, 0.906, 0.79, 0.773, 0.741, 0.656, 0.665, 0.602]'*1e11;
md.materials.density=[10925., 10925., 5506.42, 5491.45, 5456.57, 5357.06, 5307.24, 5207.13, 5156.69, 5054.69, 5002.99, 4897.83, 4844.22, 4734.6, 4678.44, 4563.07, 4503.72, 4443.16, 4412.41, 3992.14, 3983.99, 3975.84, 3912.82, 3786.78, 3723.78, 3516.39, 3489.51, 3435.78, 3359.5, 3367.1, 3184.3]';
md.materials.viscosity=[0., 0., 7.999999999999999E+21, 8.5E+21, 8.999999999999999E+21, 3.E+22, 4.E+22, 5.0000000000000004E+22, 6.E+22, 5.0000000000000004E+22, 4.5E+22, 3.E+22, 2.5000000000000002E+22, 1.7999999999999998E+22, 1.3E+22, 7.999999999999999E+21, 6.999999999999999E+21, 6.5E+21, 6.E+21, 5.5E+21, 5.E+21, 4.4999999999999995E+21, 3.9999999999999995E+21, 2.5E+21, 1.9999999999999997E+21, 1.5E+21, 9.999999999999999E+20, 6.E+20, 5.5000000000000007E+20, 2.E+20, 1.E40]';
md.materials.lame_lambda=md.materials.lame_mu*0+5e17;
md.materials.issolid=md.materials.lame_mu>0;
md.materials.numlayers=length(md.materials.lame_mu);
md.materials.burgers_mu=md.materials.lame_mu;
md.materials.burgers_viscosity=md.materials.viscosity;
md.materials.rheologymodel=md.materials.issolid*0;
md.materials.burgers_mu=md.materials.lame_mu/3;
md.materials.burgers_viscosity=md.materials.viscosity/10;
md.materials.ebm_alpha= ones(md.materials.numlayers,1)*.9;
md.materials.ebm_delta= ones(md.materials.numlayers,1)*0.2;
md.materials.ebm_taul= ones(md.materials.numlayers,1)*54*60; %54min
md.materials.ebm_tauh= ones(md.materials.numlayers,1)*18.6*cst/1e3; %18.6yr

md.love.allow_layer_deletion=1;
md.love.frequencies=[0];
md.love.nfreq=length(md.love.frequencies);
md.love.sh_nmin=1;
md.love.sh_nmax=256;
md.love.underflow_tol=1e-20;
md.love.pw_threshold=1e-3;
md.love.Gravitational_Constant=6.6732e-11;
md.love.min_integration_steps=100;
md.love.allow_layer_deletion=1;
md.love.forcing_type=11;
md.love.chandler_wobble=0;
md.love.complex_computation=0;
md.love.istemporal=1;
md.love.n_temporal_iterations=8;
%md.love.time=[0; (logspace(-3,6, 202))'*cst];
md.love.time=[0; 10.^([-3.0000e+00  -2.1045e+00  -1.2090e+00  -3.1343e-01   5.8209e-01   1.4776e+00   2.3731e+00   3.2687e+00   4.1642e+00   5.0597e+00   5.9552e+00])'*cst];
md.love=md.love.build_frequencies_from_time;

md.love.love_kernels=1;

md=solve(md,'lv');

h=md.results.LoveSolution.LoveHt;
l=md.results.LoveSolution.LoveLt;
k=md.results.LoveSolution.LoveKt;

%Fields and tolerances to track changes

field_names     ={'LoveH_loading_elastic','LoveK_loading_elastic','LoveL_loading_elastic'};
field_tolerances={2e-5,2e-5,2e-5};
field_values={h,k,l};

