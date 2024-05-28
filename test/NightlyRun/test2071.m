%Test Name: TemporalLoveNumbers
%Forward Love number solution for a viscoelastic earth, model M3-L70-V01 from
%Spada, G., Barletta, V. R., Klemann, V., Riva, R. E. M., Martinec, Z.,
%Gasperini, P., Lund, B., Wolf, D., Vermeersen, L. L. A. and King, M. A.
%(2011), A benchmark study for glacial isostatic adjustment codes. Geophysical
%Journal International, 185: 106--132. doi:10.1111/j.1365-246X.2011.04952.x

md=model();
md.cluster=generic('name',oshostname(),'np',3);

md.materials=materials('litho');
md.miscellaneous.name='test2071';
md.groundingline.migration='None';

md.verbose=verbose('all');
md.verbose=verbose('1111111111111111');
yts=365.25*24*3600;

md.materials.numlayers=6;
md.materials.radius =  [10 1222.5 3.4800e+03   5.7010e+03   5.9510e+03   6.3010e+03   6.3710e+03]'*1e3;
md.materials.density=  [1.0750e4 1.0750e+04   4.9780e+03   3.8710e+03   3.4380e+03   3.0370e+03]';
md.materials.lame_mu=  [1e-5         0   2.2834e+00   1.0549e+00   7.0363e-01   5.0605e-01]'*1e11;
md.materials.viscosity=[0            0   2.0000e+00   1.0000e+00   1.0000e+00   1.0000e+25]'*1e21;
md.materials.lame_lambda=md.materials.lame_mu*0+5e17;
md.materials.issolid=[1 0 1 1 1 1]';
md.materials.rheologymodel=zeros(md.materials.numlayers,1);
md.materials.burgers_mu=md.materials.lame_mu/3;
md.materials.burgers_viscosity=md.materials.viscosity/10;
md.materials.ebm_alpha= ones(md.materials.numlayers,1)*.9;
md.materials.ebm_delta= ones(md.materials.numlayers,1)*0.2;
md.materials.ebm_taul= ones(md.materials.numlayers,1)*54*60; %54min
md.materials.ebm_tauh= ones(md.materials.numlayers,1)*18.6*yts; %18.6yr

md.love.allow_layer_deletion=1;
md.love.frequencies=[0];
md.love.nfreq=length(md.love.frequencies);
md.love.sh_nmin=1;
md.love.sh_nmax=1000;
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
md.love.time=[0; (logspace(0,4.3, 99))'*yts];
md.love=md.love.build_frequencies_from_time;
md.love.love_kernels=1;
md.solidearth.lovenumbers.tk2secular= 0.9668;
md.solidearth.rotational.equatorialmoi=8.0131e37;
md.solidearth.rotational.polarmoi=8.0394e37;
md.solidearth.rotational.angularvelocity=7.292115e-5;

md=solve(md,'lv');

h=md.results.LoveSolution.LoveHt;
l=md.results.LoveSolution.LoveLt;
k=md.results.LoveSolution.LoveKt;

%Fields and tolerances to track changes
field_names     ={'LoveH_loading_temporal','LoveK_loading_temporal','LoveL_loading_temporal'};
field_tolerances={5e-6,5e-5,4e-5};
field_values={h,k,l};

spada=0;
if (spada)
	t=md.love.time/yts;
	addpath ../Data/
	load spada.mat
	load spada_n1.mat;
	s_weak=[9 12 15];

	hspada(:,s_weak)=[];
	kspada(:,s_weak)=[];
	lspada(:,s_weak)=[];

	hts=zeros(length(t),257);
	lts=zeros(length(t),257);
	kts=zeros(length(t),257);

	d=1;
	hts(:,d+1)=hspada_n1(2);
	lts(:,d+1)=lspada_n1(2);
	kts(:,d+1)=-1;
	for mo=1:9 
		hts(:,d+1)=hts(:,d+1)-hspada_n1(3+mo)./sspada_n1(1+mo).*(1-exp(t/1e3*sspada_n1(1+mo)));
		lts(:,d+1)=lts(:,d+1)-lspada_n1(3+mo)./sspada_n1(1+mo).*(1-exp(t/1e3*sspada_n1(1+mo)));
	end

	for d=2:256
		hts(:,d+1)=hspada(d-1,2);
		lts(:,d+1)=lspada(d-1,2);
		kts(:,d+1)=kspada(d-1,2);
		for mo=1:9 
			hts(:,d+1)=hts(:,d+1)-hspada(d-1,3+mo)./sspada(d-1,1+mo).*(1-exp(t/1e3*sspada(d-1,1+mo)));
			lts(:,d+1)=lts(:,d+1)-lspada(d-1,3+mo)./sspada(d-1,1+mo).*(1-exp(t/1e3*sspada(d-1,1+mo)));
			kts(:,d+1)=kts(:,d+1)-kspada(d-1,3+mo)./sspada(d-1,1+mo).*(1-exp(t/1e3*sspada(d-1,1+mo)));
		end
	end
end

return


h=md.results.LoveSolution.LoveHt;
l=md.results.LoveSolution.LoveLt;
k=md.results.LoveSolution.LoveKt;
ht=h';
lt=l';
kt=k';
th=md.results.LoveSolution.LoveTidalHt';
tl=md.results.LoveSolution.LoveTidalLt';
tk=md.results.LoveSolution.LoveTidalKt';
tht=ht*0;tht(3,:)=th(3,:);
tkt=kt*0;tkt(3,:)=tk(3,:);
tlt=lt*0;tlt(3,:)=tl(3,:);
pmtf1=md.results.LoveSolution.LovePMTF1t(:,3)';
pmtf2=md.results.LoveSolution.LovePMTF2t(:,3)';
time=md.love.time/yts;

save ../Data/lnb_temporal ht kt lt tht tkt tlt pmtf1 pmtf2 time;

