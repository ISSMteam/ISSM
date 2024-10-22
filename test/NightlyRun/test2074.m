%Test Name: EBMLoveNumbers

md=model();
md.cluster=generic('name',oshostname(),'np',3);

md.materials=materials('litho');
md.miscellaneous.name='test2074';

md.verbose=verbose('all');
md.verbose=verbose('1111111111111111');
yts=365.25*24*3600;

md.materials.numlayers=6;
md.materials.radius =  [100 1222.5 3.4800e+03   5.7010e+03   5.9510e+03   6.3310e+03   6.3710e+03]'*1e3;
md.materials.density=  [1.0750e4 1.0750e+04   4.9780e+03   3.8710e+03   3.4380e+03   3.0370e+03]';
md.materials.lame_mu=  [1e-5         0   2.2834e+00   1.0549e+00   7.0363e-01   5.0605e-01]'*1e11;
md.materials.viscosity=[0            0   20.0000e+00   0.5000e+00   0.078600e+00   1.0000e+25]'*1e21;
md.materials.lame_lambda=md.materials.lame_mu*0+5e17;
md.materials.issolid=[1 0 1 1 1 1]';
md.materials.rheologymodel=zeros(md.materials.numlayers,1)+2;
md.materials.rheologymodel(end)=0; %let the lithosphere be maxwell so we dont end up with a small unrelaxed mu for large timescales
md.materials.ebm_alpha= ones(md.materials.numlayers,1)*.4;
md.materials.ebm_delta= ones(md.materials.numlayers,1)*.9;
md.materials.ebm_taul= ones(md.materials.numlayers,1)*54*60; %54min
md.materials.ebm_tauh= ones(md.materials.numlayers,1)*7.134*yts; %7.134yr
md.materials.setlitho2prem;
md.materials.lame_lambda=md.materials.lame_mu*0+5e17;

md.love.frequencies=[0];
md.love.nfreq=length(md.love.frequencies);
md.love.sh_nmin=1;
md.love.sh_nmax=256;
md.love.underflow_tol=1e-20;
md.love.Gravitational_Constant=6.6732e-11;
md.love.min_integration_steps=100;
md.love.max_integration_dr=5e3;

md.love.istemporal=1;
md.love.n_temporal_iterations=6;
md.love.time=[0; (logspace(0,4.3, 24))'*yts];

load ../Data/hypergeom.mat;
md.love.hypergeom_table1=h1real;
md.love.hypergeom_table2=h1real*0;
md.love.hypergeom_nalpha=101;
md.love.hypergeom_nz=length(z);
md.love.hypergeom_z=z;

if md.love.istemporal
	md.love=md.love.build_frequencies_from_time;
end

md.solidearth.lovenumbers.tk2secular= 0.9668;
md.solidearth.rotational.equatorialmoi=8.0131e37;
md.solidearth.rotational.polarmoi=8.0394e37;
md.solidearth.rotational.angularvelocity=7.292115e-5;

md=solve(md,'lv');
hebm=md.results.LoveSolution.LoveHt;
kebm=md.results.LoveSolution.LoveKt;
lebm=md.results.LoveSolution.LoveLt;

%Fields and tolerances to track changes

field_names     ={'LoveHebm','LoveKebm','LoveLebm'};
field_tolerances={1.0e-6,1.0e-6,1.0e-6};
field_values={hebm,kebm,lebm};