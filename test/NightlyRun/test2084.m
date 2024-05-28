%Test Name: GiaCaron
%Forward Love number solution for a viscoelastic earth, model M3-L70-V01 from
%Spada, G., Barletta, V. R., Klemann, V., Riva, R. E. M., Martinec, Z.,
%Gasperini, P., Lund, B., Wolf, D., Vermeersen, L. L. A. and King, M. A.
%(2011), A benchmark study for glacial isostatic adjustment codes. Geophysical
%Journal International, 185: 106--132. doi:10.1111/j.1365-246X.2011.04952.x

md=model();
md.cluster=generic('name',oshostname(),'np',8);

% set validation=1 for comparing against the Spada benchmark.
validation=0; 

md.materials=materials('litho');
md.miscellaneous.name='FourierLoveTest';
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
md.materials.burgers_mu=md.materials.lame_mu/3;
md.materials.burgers_viscosity=md.materials.viscosity/10;
md.materials.ebm_alpha= ones(md.materials.numlayers,1)*.9;
md.materials.ebm_delta= ones(md.materials.numlayers,1)*0.2;
md.materials.ebm_taul= ones(md.materials.numlayers,1)*54*60; %54min
md.materials.ebm_tauh= ones(md.materials.numlayers,1)*18.6*cst/1e3; %18.6yr
%setlitho2prem(md.materials)

md.love.allow_layer_deletion=1;
md.love.frequencies=[0; (logspace(-6,3, 1000))'/cst];
md.love.nfreq=length(md.love.frequencies);
md.love.sh_nmin=1;
md.love.sh_nmax=1000;
md.love.underflow_tol=1e-20;
md.love.pw_threshold=1e-3;
md.love.Gravitational_Constant=6.6732e-11;
md.love.allow_layer_deletion=1;
md.love.forcing_type=11;
md.love.chandler_wobble=0;
md.love.complex_computation=0;

md.love.istemporal=1;
md.love.n_temporal_iterations=8;
%md.love.time=(logspace(-6,5, 2))'*cst;
md.love.time=[0; (logspace(-3,5, 24))'*cst];

%md.love.time=(linspace(1/12,10, 10*12))'*cst/1e3;
md.love.love_kernels=1;
if md.love.istemporal
	md.love=md.love.build_frequencies_from_time;
end

md=solve(md,'lv');

ht=md.results.LoveSolution.LoveHt';
lt=md.results.LoveSolution.LoveLt';
kt=md.results.LoveSolution.LoveKt';
t=md.love.time/cst*1e3;

%hs=reshape(md.results.LoveSolution.LoveHi(:,:), [ md.love.sh_nmax+1, 2*md.love.n_temporal_iterations, length(t)]);
%hs=permute(hs,[3 2 1]);
%[ht,h,hsig,hconv]=postwidder_love(md,md.love.n_temporal_iterations,t,hs,1e-5);

%Fields and tolerances to track changes
%loading love numbers
field_names     ={'LoveH_loading_elastic','LoveK_loading_elastic','LoveL_loading_elastic'};
field_tolerances={2.0e-8,2.0e-8,2.0e-8};
field_values={...
	(md.results.LoveSolution.LoveHt(:,1)),...
	(md.results.LoveSolution.LoveKt(:,1)),...
	(md.results.LoveSolution.LoveLt(:,1)),...
	};

return

if false
addpath('../../../../invGIA/Spada_benchmark/')
load spada.mat
s_weak=[9 12 15];

hspada(:,s_weak)=[];
kspada(:,s_weak)=[];
lspada(:,s_weak)=[];

hts=zeros(length(t),md.love.sh_nmax+1);
lts=zeros(length(t),md.love.sh_nmax+1);
kts=zeros(length(t),md.love.sh_nmax+1);
for d=max(2,md.love.sh_nmin):md.love.sh_nmax
hts(:,d+1)=hspada(d-1,2);
lts(:,d+1)=lspada(d-1,2);
kts(:,d+1)=kspada(d-1,2);
for mo=1:9 
hts(:,d+1)=hts(:,d+1)-hspada(d-1,3+mo)./sspada(d-1,1+mo).*(1-exp(t/1e3*sspada(d-1,1+mo)));
lts(:,d+1)=lts(:,d+1)-lspada(d-1,3+mo)./sspada(d-1,1+mo).*(1-exp(t/1e3*sspada(d-1,1+mo)));
kts(:,d+1)=kts(:,d+1)-kspada(d-1,3+mo)./sspada(d-1,1+mo).*(1-exp(t/1e3*sspada(d-1,1+mo)));
end
end
close all

if md.love.sh_nmin==md.love.sh_nmax
subplot(2,3,1)
plot(t/cst*1e3,ht(:,md.love.sh_nmin+1),'x-',t/cst*1e3,hts);
set(gca, 'xscale', 'log'); shading flat;
title('h')

subplot(2,3,2)
plot(t/cst*1e3,lt(:,md.love.sh_nmin+1),'x-',t/cst*1e3,lts);
set(gca, 'xscale', 'log'); shading flat;
title('l')

subplot(2,3,3)
plot(t/cst*1e3,kt(:,md.love.sh_nmin+1),'x-',t/cst*1e3,kts);
set(gca, 'xscale', 'log'); shading flat;
title('k')

subplot(2,3,4)
plot(t/cst*1e3,log10(abs((ht-hts)./hts))');
set(gca, 'xscale', 'log'); shading flat;

subplot(2,3,5)
plot(t/cst*1e3,log10(abs((lt-lts)./lts))');
set(gca, 'xscale', 'log'); shading flat;
title('l')

subplot(2,3,6)
plot(t/cst*1e3,log10(abs((kt-kts)./kts))');
set(gca, 'xscale', 'log'); shading flat;
title('k')
else

[T,N]=meshgrid(t/cst*1e3,0:md.love.sh_nmax);
subplot(1,3,1)
pcolor(T,N,log10(abs((ht-hts)./hts))');
set(gca, 'xscale', 'log'); shading flat;colorbar;
title('h')

subplot(1,3,2)
pcolor(T,N,log10(abs((lt-lts)./lts))');
set(gca, 'xscale', 'log'); shading flat;colorbar;
title('l')

subplot(1,3,3)
pcolor(T,N,log10(abs((kt-kts)./kts))');
set(gca, 'xscale', 'log'); shading flat;colorbar;
title('k')
end


% validate elastic loading solutions against the Spada benchmark. {{{ 
if validation 
	spada_solutions = load('spada_elastic_loading_deg_h_l_k'); 
	spada_d = spada_solutions(:,1); 
	spada_h = spada_solutions(:,2); 
	spada_l = spada_solutions(:,3); 
	spada_k = spada_solutions(:,4); 

	%rename ISSM solutions.  
	issm_d = [1:md.love.sh_nmax];  
	issm_h = md.results.LoveSolution.LoveHr(2:end,1); 
	issm_l = md.results.LoveSolution.LoveLr(2:end,1); 
	issm_k = md.results.LoveSolution.LoveKr(2:end,1); 

	% relative difference for each degree, except for zero.  
	diff_h = 1 - issm_h./spada_h;  
	diff_l = 1 - issm_l./spada_l;  
	diff_k = 1 - issm_k./spada_k;  

	figure 
	plot(spada_d,[diff_h diff_l diff_k]); grid on; 
	legend('h','l','k'); title('loading love'); 

else
	% 
end 

end
% }}} 
