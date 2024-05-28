%Test Name: SESAWslc
% SESAW method of solving GRD slc
% reference: Adhikari et al., 2016, GMD, https://doi.org/10.5194/gmd-9-1087-2016

md=model;
load ../Data/SlcTestMesh.mat;
md.mesh=SlcMesh; %700 km resolution mesh

% read in love numbers.
love_numbers = lovenumbers('maxdeg',10000);

% compute Green's functions
disp(['Computing Greens functions...']);
[Grigid,Gelastic,Uelastic]=greensfunctions(md.mesh.elements,md.mesh.lat,md.mesh.long,love_numbers);
greens.Grigid = Grigid;
greens.Gelast = Gelastic;
greens.Uelast = Uelastic;
clearvars Grigid Gelastic Uelastic

% compute lat,long at elemental centroids.
[late,longe] = latelonge(md.mesh.elements,md.mesh.lat,md.mesh.long);

% load GRACE data. Same as used in test2020.
load('../Data/GRACE_JPL_April2002_WEH.mat');
lat = repmat(lat',720,1);
lon = repmat(lon,1,360);
F = scatteredInterpolant(lat(:),lon(:),weh(:));

% map GRACE data onto elemental centorids.
loads_element = F(late,longe);
loads_element(isnan(loads_element))=0;

% ocean mask mapped onto the elemental centroids.
ocean_element = gmtmask(late,longe);

% Area of individual elements
area_element=GetAreasSphericalTria(md.mesh.elements,md.mesh.lat,md.mesh.long,md.solidearth.planetradius);

% Parameters input for SESAWslc solver.
para.ocean_element = ocean_element;
para.loads_element = loads_element;
para.area_element = area_element;
para.earth_density = md.materials.earth_density;
para.ocean_density = md.materials.rho_water;
para.loads_density = md.materials.rho_freshwater; % if land loads are ice, use ice density.

para.rel_tol = 1e-5;

% solid earth rheology.
para.solidearth = 'rigid'; % 'rigid' or 'elastic';

% rotational feedbacks.
para.rotational.flag = 0; % Rotational flag on (1) or off (0)
para.rotational.earth_radius = md.solidearth.planetradius;
para.rotational.load_love_k2 = love_numbers.k(3);
para.rotational.tide_love_k2 = love_numbers.tk(3);
para.rotational.tide_love_h2 = love_numbers.th(3);
para.rotational.tide_love_k2secular = love_numbers.tk2secular;
para.rotational.moi_p = md.solidearth.rotational.polarmoi;
para.rotational.moi_e = md.solidearth.rotational.equatorialmoi;
para.rotational.omega = md.solidearth.rotational.angularvelocity;

% solve: Rigid without rotational feedbacks.
disp(['Solving sesaw-slc for Rigid Earth WITHOUT rotational feedback...']);
[eus_rigid,rsl_rigid] = SESAWslr(md.mesh.elements,md.mesh.lat,md.mesh.long,greens,para);

% solve: Rigid with rotational feedbacks.
para.rotational.flag = 1;
disp(['Solving sesaw-slc for Rigid Earth WITH rotational feedback...']);
[eus_rigid_rot,rsl_rigid_rot] = SESAWslr(md.mesh.elements,md.mesh.lat,md.mesh.long,greens,para);

% solve: Elastic with rotational feedbacks.
para.solidearth = 'elastic';
disp(['Solving sesaw-slc for Elastic Earth WITH rotational feedback...']);
[eus_elast_rot,rsl_elast_rot] = SESAWslr(md.mesh.elements,md.mesh.lat,md.mesh.long,greens,para);

% solve: Elastic with rotational feedbacks.
para.rotational.flag = 0;
disp(['Solving sesaw-slc for Elastic Earth WITHOUT rotational feedback...']);
[eus_elast,rsl_elast] = SESAWslr(md.mesh.elements,md.mesh.lat,md.mesh.long,greens,para);

%Fields and tolerances to track changes
field_names={'eus_rigid','eus_rigid_rot','eus_elast','eus_elast_rot',...
             'rsl_rigid','rsl_rigid_rot','rsl_elast','rsl_elast_rot'};
field_tolerances={1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13};
field_values={eus_rigid,eus_rigid_rot,eus_elast,eus_elast_rot,...
             rsl_rigid,rsl_rigid_rot,rsl_elast,rsl_elast_rot};

