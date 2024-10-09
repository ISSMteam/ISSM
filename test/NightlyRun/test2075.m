%Test Name: Complex_Lovenumbers
%Forward Love number solution for a viscoelastic earth, model M3-L70-V01 from
%Spada, G., Barletta, V. R., Klemann, V., Riva, R. E. M., Martinec, Z.,
%Gasperini, P., Lund, B., Wolf, D., Vermeersen, L. L. A. and King, M. A.
%(2011), A benchmark study for glacial isostatic adjustment codes. Geophysical
%Journal International, 185: 106--132. doi:10.1111/j.1365-246X.2011.04952.x

md=model();
md.cluster=generic('name',oshostname(),'np',1);

md.materials=materials('litho');
md.miscellaneous.name='test2075';
md.groundingline.migration='None';

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
md.materials.ebm_taul= ones(md.materials.numlayers,1)*54*60*100; %54min
md.materials.ebm_tauh= ones(md.materials.numlayers,1)*7.134*yts; %7.134yr
md.materials.setlitho2prem;
md.materials.lame_lambda=md.materials.lame_mu*0+5e17;

md.love.frequencies=[0 logspace(-8,3,100)/yts];
md.love.nfreq=length(md.love.frequencies);
md.love.sh_nmin=2;
md.love.sh_nmax=2;
md.love.underflow_tol=1e-20;
md.love.pw_threshold=1e-3;
md.love.Gravitational_Constant=6.6732e-11;
md.love.min_integration_steps=500;
md.love.max_integration_dr=5e3;
md.love.complex_computation=1;

md.love.istemporal=0;
md.love.time=[];


load ../Data/hypergeom.mat;
md.love.hypergeom_table1=h1complex;
md.love.hypergeom_table2=h2complex;
md.love.hypergeom_nalpha=101;
md.love.hypergeom_nz=length(z);
md.love.hypergeom_z=z;

if md.love.istemporal
	md.love=md.love.build_frequencies_from_time;
end
md.love.love_kernels=0;
md.solidearth.lovenumbers.tk2secular= 0.9668;
md.solidearth.rotational.equatorialmoi=8.0131e37;
md.solidearth.rotational.polarmoi=8.0394e37;
md.solidearth.rotational.angularvelocity=7.292115e-5;

md=solve(md,'lv');
krebm=1+md.results.LoveSolution.LoveKf(:,3);
kiebm=1+md.results.LoveSolution.LoveKfi(:,3)-1;

md.materials.rheologymodel=zeros(md.materials.numlayers,1);
md=solve(md,'lv');
krmax=1+md.results.LoveSolution.LoveKf(:,3);
kimax=1+md.results.LoveSolution.LoveKfi(:,3)-1;

clf
f=md.love.frequencies*yts;
semilogx(f,krmax,f,krebm);
hold on
semilogx(f,kimax,'--',f,kiebm,'--');

return
subplot(1,2,1)
plot(1:90, kebm(tt,2:end)./kmax(1,2:end));
%ylim([1 2.5])
xlabel('degree')
set(gca, 'Fontsize', 18')
ylabel('$$\frac{1+k(t)_{EBM}}{1+k_e}$$', 'interpreter', 'LaTeX')
legend('t=1 month', 't=1 year', 't=5 years', 't=20 years')
subplot(1,2,2)                           
plot(1:90, kebm(tt,2:end)./kmax(tt,2:end))
%ylim([1 2.5])
xlabel('degree')
set(gca, 'Fontsize', 18')
ylabel('$$\frac{1+k(t)_{EBM}}{1+k(t)_{Maxwell}}$$', 'interpreter', 'LaTeX')
legend('t=1 month', 't=1 year', 't=5 years', 't=20 years')

figure

plot(md.love.time(2:end-1)/yts, (kebm(3:end,3)-kebm(1:end-2,3))./(md.love.time(3:end)-md.love.time(1:end-2))*yts/(kmax(1,3)));
ylabel('$$\frac{\dot{k_2}(t)_{EBM}}{1+k_2(t)_{Elastic}}$$', 'interpreter', 'LaTeX')
set(gca, 'Fontsize', 18')
xlabel('t (yr)');
