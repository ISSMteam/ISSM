%Test Name: ThermalConduction
% This file can be run to check that the conduction is correctly modeled.
% There is no velocity (no advection) the only thermal boundary conditions are an imposed temperature
% at the lower and upper surface. The result must be a linear temperature from the upper to the lower
% surface. if it is not the case, something is thermal modeling has been changed...
printingflag=false;

md=model();
md=triangle(md,'../Exp/Square.exp',100000.);
md=setmask(md,'all','');
md=parameterize(md,'../Par/SquareThermal.par');
md=extrude(md,11,2.);
md=setflowequation(md,'HO','all');

pos1=find(md.mesh.elementonbase);     md.thermal.spctemperature(md.mesh.elements(pos1,1:3))=10.;
pos2=find(md.mesh.elementonsurface); md.thermal.spctemperature(md.mesh.elements(pos2,4:6))=0.;
md.initialization.pressure=zeros(md.mesh.numberofvertices,1);

%analytical results
%d2T/dz2=0 T(bed)=10 T(surface)=0  => T=0*(z-bed)/thickness+10*(surface-z)/thickness
%each layer of the 3d mesh must have a constant value
md.initialization.temperature=10.*(md.geometry.surface-md.mesh.z)./md.geometry.thickness;

%modeled results
md.cluster=generic('name',oshostname(),'np',2);
md=solve(md,'Thermal');

%plot results
comp_temp=md.results.ThermalSolution.Temperature;
relative=abs((comp_temp-md.initialization.temperature)./md.initialization.temperature)*100.;
relative(find(comp_temp==md.initialization.temperature))=0.;
plotmodel(md,'data',comp_temp,'title','Modeled temperature [K]','data',md.initialization.temperature,'view',3,...
	'title','Analytical temperature [K]','view',3,'data',comp_temp-md.initialization.temperature,...
	'title','Absolute error [K]','view',3,'data',relative,'title','Relative error [%]','view',3,...
	'figposition','mathieu','FontSize#all',20)
if printingflag,
	set(gcf,'Color','w')
	printmodel('thermalconduction','png','margin','on','marginsize',25,'frame','off','resolution',0.7,'hardcopy','off');
	system(['mv thermalconduction.png ' ISSM_DIR '/website/doc_pdf/validation/Images/Thermal ']);
end

%Fields and tolerances to track changes
field_names     ={'ConductionTemperature'};
field_tolerances={1e-13};
field_values    ={comp_temp};
