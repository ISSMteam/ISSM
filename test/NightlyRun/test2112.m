%Test Name: Esa2Dsurface
%AIS 2 -- load centered at the south pole! 

%mesh ais: {{{
md=model();
md=triangle(md,'../Exp/Ais.exp',100000); % max element size
% }}}
%define load: {{{
md.esa.deltathickness=zeros(md.mesh.numberofelements,1);
disc_radius=500; % km
index=md.mesh.elements;
x_element=mean(md.mesh.x(index),2); 
y_element=mean(md.mesh.y(index),2); 
rad_dist=sqrt(x_element.^2+y_element.^2)/1000;  % radial distance in km
pos=find(rad_dist<=disc_radius);
md.esa.deltathickness(pos)=-1;   % 1 m water withdrawal
% }}}
%read in love numbers:{{{
md.solidearth.lovenumbers=lovenumbers('maxdeg',10000,'referenceframe','CF');
% }}}
%mask:  {{{
%make sure wherever there is an ice load, that the mask is set to ice: 
md.mask.ice_levelset=ones(md.mesh.numberofvertices,1);
pos=find(md.esa.deltathickness); md.mask.ice_levelset(md.mesh.elements(pos,:))=-1;

%is ice grounded? 
md.mask.ocean_levelset=-ones(md.mesh.numberofvertices,1);
pos=find(md.mask.ice_levelset<=0); md.mask.ocean_levelset(pos)=1; 
% }}}
%geometry:  {{{
di=md.materials.rho_ice/md.materials.rho_water;
md.geometry.thickness=ones(md.mesh.numberofvertices,1);
md.geometry.surface=(1-di)*zeros(md.mesh.numberofvertices,1);
md.geometry.base=md.geometry.surface-md.geometry.thickness;
md.geometry.bed=md.geometry.base;
% }}}
%materials:  {{{
md.initialization.temperature=273.25*ones(md.mesh.numberofvertices,1);
md.materials.rheology_B=paterson(md.initialization.temperature);
md.materials.rheology_n=3*ones(md.mesh.numberofelements,1);
% }}}
%additional parameters, miscellaneous: {{{
md.miscellaneous.name='test2112';
md.esa.degacc=0.01;
md.esa.hemisphere = -1; 
% }}}

%solve esa: {{{
md.esa.requested_outputs = {'EsaUmotion','EsaNmotion','EsaEmotion','EsaXmotion','EsaYmotion'}; 
md.cluster=generic('name',oshostname(),'np',3); 
md.verbose=verbose('111111111');
md=solve(md,'Esa');
% }}}
%fields and tolerances to track changes {{{
field_names     ={'EsaUmotion','EsaNmotion','EsaEmotion','EsaXmotion','EsaYmotion'}; 
field_tolerances={1e-13,4e-13,3e-12,3e-13,3e-13};
field_values={...
	(md.results.EsaSolution.EsaUmotion),...
	(md.results.EsaSolution.EsaNmotion),...
	(md.results.EsaSolution.EsaEmotion),...
	(md.results.EsaSolution.EsaXmotion),...
	(md.results.EsaSolution.EsaYmotion),...
	};
% }}} 

