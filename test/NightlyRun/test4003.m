%Test Name: IceOceanCoupling Dan Goldberg'd setup
%ISSM/MITgcm coupled set-up
%
%Script control parameters
steps=1:11;

%To download and recompile MITgcm from scratch:
!rm -rf ${ISSM_DIR}/test/MITgcm/install
!rm -rf ${ISSM_DIR}/test/MITgcm/build/*
!rm -rf Models

%Organizer
!mkdir Models
org=organizer('repository','Models','prefix','IceOcean.','steps',steps);

presentdirectory=pwd;

% {{{ Parameters:
if perform(org,'Parameters'),
    Nx=3;   % number of longitude cells
    Ny=200; % number of latitude cells
    Nz=90;  % number of MITgcm vertical cells
    nPx=1;  % number of MITgcm processes to use in x direction
    nPy=8;  % number of MITgcm processes to use in y direction
    xgOrigin=0;     % origin of longitude
    ygOrigin=-75.5; % origin of latitude
    dLong=.125;     % longitude grid spacing
    dLat=dLong/16;  % latitude grid spacing
    delZ=10;        % thickness of vertical levels (m)
    gravity= 9.81;  % gravity (m^2/s)
    rho_ice=917;
    rho_water=1030;
    di=rho_ice/rho_water;
    prec = 'real*8'; % precision of MITgcm input binary files

    % bathymetry and ice sheet geometry
    H = -900;	    % water depth in the ice shelf cavity
    Hmin = -600;    % deepest point of cavern		
    Hmax = -300;    % shallowest point of cavern
    jEnd = Ny*3/4;  % where ice-shelf ends
    j2   = jEnd+1;
    
    % initial ocean conditions
    T_sfc = -2;
    T_bot = -1.9;
    S_sfc = 34.2;
    S_bot = 34.3;

    savedata(org, Nx, Ny, nPx, nPy, Nz, dLong, dLat, delZ, xgOrigin, ygOrigin, ...
             rho_ice, rho_water, di, H, Hmin, Hmax, jEnd, j2, gravity, prec, ...
             T_sfc,T_bot,S_sfc,S_bot);
end
% }}}
% {{{ Bathymetry:
if perform(org,'Bathymetry'),

    loaddata(org,'Parameters');

    %create lat,lon
    latg = ygOrigin+[0:Ny-1]*dLat;
    latc = latg+dLat/2;
    long = xgOrigin+[0:Nx-1]*dLong;
    lonc = long+dLong/2;
    [lat lon]=meshgrid(latc,lonc);
    zC=-delZ*([1:Nz]-0.5);
    zF=-delZ*[0:Nz];
    
    %create bathymetry:
    bathymetry = ones(Nx,Ny)*H;
    bathymetry(:,end) = 0;

    %save bathymetry file for MITgcm
    savedata(org,lat,lon,bathymetry);

end
% }}}
% {{{ IceSheetGeometry:
if perform(org,'IceSheetGeometry'),

    loaddata(org,'Parameters');
    loaddata(org,'Bathymetry');
    latmin=min(lat(:));
    latmax=max(lat(:));

    dHdy = (Hmax-Hmin)/dLat/(jEnd-2); %Slope of ice shelf
    draft=bathymetry;
    for i=1:Nx
        draft(i,:)=Hmin+dHdy*[-1:Ny-2]*dLat;
    end
    draft(:,j2:Ny)=0;

    ice_mask=ones(Nx,Ny);
    ice_mask(:,j2:Ny)=0;
    iceshelf_mask=ice_mask;
    thickness=abs(draft)/di;
    
    savedata(org,ice_mask,iceshelf_mask,draft,thickness);
    
    close all, figure(2), clf
    subplot(411), pcolorcen(bathymetry); colorbar, title('bathymetry')
    subplot(412), pcolorcen(ice_mask); colorbar, title('ice and iceshelf mask')
    subplot(413), pcolorcen(draft); colorbar, title('draft')
    subplot(414), pcolorcen(thickness); colorbar, title('thickness')
    
end
% }}}

%Configure MITgcm
% {{{ GetMITgcm:
if perform(org,'GetMITgcm'),
  system([pwd '/../MITgcm/get_mitgcm.sh']);
end
% }}}
% {{{ BuildMITgcm:
if perform(org,'BuildMITgcm'),
    system(['../MITgcm/build_4003.sh generic ' pwd '/../MITgcm']);
end
% }}}
addpath(recursivepath([pwd '/../MITgcm']));
% {{{ RunUncoupledMITgcm:
if perform(org,'RunUncoupledMITgcm'),

    loaddata(org,'Parameters');
    loaddata(org,'Bathymetry');
    loaddata(org,'IceSheetGeometry');

    % rename previous run directory and create new one
    if exist ('run.old')
        !\rm -rf run.old
    end
    if exist ('run')
        !\mv run run.old
    end
    !\mkdir run
    !\cp ../MITgcm/build/mitgcmuv run
    !\cp ../MITgcm/input_4003/* run
    !\cp ../MITgcm/input_4003/eedata_uncoupled run/eedata
    
    % {{{ Construct MITgcm binary input files
    namF='run/bathy_flat.bin';
    fid=fopen(namF,'w','b'); fwrite(fid,bathymetry,prec);fclose(fid);
    
    namF='run/shelficeTopo.Lin.bin';
    fid=fopen(namF,'w','b'); fwrite(fid,draft,prec);fclose(fid);
    
    del_T = (T_bot - T_sfc)/(59*delZ);
    tref=zeros(1,Nz);
    for k = 1:Nz;
        tref(k) = T_sfc + del_T*((k-20)*delZ);
        tref(k)= max(T_sfc,min(tref(k),T_bot));
    end
    namF='run/temp_obc.bin';
    tref=[tref; tref; tref];
    fid=fopen(namF,'w','b'); fwrite(fid,tref,prec);fclose(fid);

    del_S = (S_bot - S_sfc)/(59*delZ);
    sref=zeros(1,Nz);
    for k = 1:Nz;
        sref(k) = S_sfc + del_S*((k-20)*delZ);
        sref(k)= max(S_sfc,min(sref(k),S_bot));
    end
    namF='run/salt_obc.bin';
    sref=[sref; sref; sref];
    fid=fopen(namF,'w','b'); fwrite(fid,sref,prec);fclose(fid);

    zax=[1:Nz];
    v1=2.5e-2;
    var=1+Nz-2*zax; var=var/(Nz-1);
    vobc=v1*var;
    namF='run/vVel_obc.bin';
    vobc=[vobc; vobc; vobc];
    fid=fopen(namF,'w','b'); fwrite(fid,vobc,prec);fclose(fid);

    var=zeros(Nx,Ny,Nz);
    for i=1:Nx, for j=1:Ny
            var(i,j,:)=tref(1,:);
        end, end
    namF='run/temp_ini.bin';
    fid=fopen(namF,'w','b'); fwrite(fid,var,prec);fclose(fid);

    for i=1:Nx, for j=1:Ny
            var(i,j,:)=sref(1,:);
        end, end
    namF='run/salt_ini.bin';
    fid=fopen(namF,'w','b'); fwrite(fid,var,prec);fclose(fid);
    % }}}

    cd run
    eval(['!mpirun -np ' int2str(nPx*nPy) ' ./mitgcmuv']);
    cd ..
end
% }}}

%Configure ISSM
% {{{ CreateMesh:
if perform(org,'CreateMesh'),

    loaddata(org,'Parameters');

    %create model:
    md=model();

    %Grab lat,long from MITgcm:
    long=readbin('run/XG.data',[Nx Ny]);
    long=[long long(:,end)]; long=[long; 2*long(Nx,:)-long(Nx-1,:)];
    lat=readbin('run/YG.data',[Nx Ny]);
    lat=[lat 2*lat(:,Ny)-lat(:,Ny-1)]; lat=[lat; lat(end,:)];

    %project lat,long:
    [x,y]=ll2xy(lat(:),long(:),-1);
    index=[];
    %  C  D
    %  A  B
    for j=1:Ny,
        for i=1:Nx,
            A=(j-1)*(Nx+1)+i;
            B=(j-1)*(Nx+1)+i+1;
            C=j*(Nx+1)+i;
            D=j*(Nx+1)+i+1;
            index(end+1,:)=[A B C];
            index(end+1,:)=[C B D];
        end
    end

    %fill mesh and model:
    md=meshconvert(md,index,x,y);
    md.mesh.lat=lat(:);
    md.mesh.long=long(:);

    savemodel(org,md);

end
% }}}
% {{{ MeshGeometry:
if perform(org,'MeshGeometry'),

    loaddata(org,'Parameters');
    loaddata(org,'CreateMesh');
    loaddata(org,'Bathymetry');
    loaddata(org,'IceSheetGeometry');

    %transfer to vertices:
    bathymetry=[bathymetry bathymetry(:,end)]; bathymetry=[bathymetry(1,:); bathymetry];
    ice_mask=[ice_mask ice_mask(:,end)]; ice_mask=[ice_mask(1,:); ice_mask];
    iceshelf_mask=[iceshelf_mask iceshelf_mask(:,end)]; iceshelf_mask=[iceshelf_mask(1,:); iceshelf_mask];
    thickness=[thickness thickness(:,end)]; thickness=[thickness; thickness(end,:)];

    %start filling some of the fields
    md.geometry.bed=bathymetry(:);
    md.geometry.thickness=thickness(:);
    md.geometry.base=-917/1028*md.geometry.thickness;
    md.geometry.surface=md.geometry.base+md.geometry.thickness;

    %nothing passes icefront:
    pos=find((~ice_mask(:) & ice_mask(:)~=0) | thickness(:)==0);
    md.geometry.thickness(pos)=1;
    md.geometry.surface(pos)=(1-di)*md.geometry.thickness(pos);
    md.geometry.base(pos)=-di*md.geometry.thickness(pos);

    %level sets:
    md.mask.ice_levelset=iceshelf_mask(:);
    pos=find(md.mask.ice_levelset==1); md.mask.ice_levelset(pos)=-1;
    pos=find(md.mask.ice_levelset==0); md.mask.ice_levelset(pos)=1;
    md.mask.ocean_levelset=-ones(md.mesh.numberofvertices,1);

    savemodel(org,md);

end
% }}}
% {{{ ParameterizeIce:
if perform(org,'ParameterizeIce'),

	loaddata(org,'MeshGeometry');

	%miscellaneous
	md.miscellaneous.name='test4003';

	%initial velocity:
	md.initialization.vx=zeros(md.mesh.numberofvertices,1);
	md.initialization.vy=zeros(md.mesh.numberofvertices,1);
	md.initialization.vz=zeros(md.mesh.numberofvertices,1);

	%friction:
	md.friction.coefficient=0*ones(md.mesh.numberofvertices,1);
	pos=find(md.mask.ocean_levelset>0);
	md.friction.coefficient(pos)=5;
	md.friction.p=ones(md.mesh.numberofelements,1);
	md.friction.q=ones(md.mesh.numberofelements,1);

	%temperatures and surface mass balance:
	md.initialization.temperature=(273.15-22)*ones(md.mesh.numberofvertices,1);
	md.initialization.pressure=md.materials.rho_ice*md.constants.g*(md.geometry.surface-md.geometry.base);
	md.smb.mass_balance = 0*ones(md.mesh.numberofvertices,1);

	%Flow law
	md.materials.rheology_B=paterson(md.initialization.temperature);
	md.materials.rheology_n=3*ones(md.mesh.numberofelements,1);
	md.damage.D=zeros(md.mesh.numberofvertices,1);
	md.damage.spcdamage=NaN*ones(md.mesh.numberofvertices,1);

	%the spcs going
	md.stressbalance.spcvx=NaN*ones(md.mesh.numberofvertices,1);
	md.stressbalance.spcvy=NaN*ones(md.mesh.numberofvertices,1);
	md.stressbalance.spcvz=NaN*ones(md.mesh.numberofvertices,1);
	md.stressbalance.referential=NaN*ones(md.mesh.numberofvertices,6);
	md.stressbalance.loadingforce=0*ones(md.mesh.numberofvertices,3);
	md.masstransport.spcthickness=NaN*ones(md.mesh.numberofvertices,1);

	%get some flux at the ice divide:
	pos=find(md.mesh.lat==min(md.mesh.lat));
	md.masstransport.spcthickness(pos)=md.geometry.thickness(pos);
	md.stressbalance.spcvy(pos)=800;
	md.stressbalance.spcvx(pos)=0;

	%deal with boundaries, excluding icefront:
	pos=find(md.mesh.long==min(md.mesh.long) | md.mesh.long==max(md.mesh.long));
	md.stressbalance.spcvx(pos)=0;

	point1=find(md.mesh.y==min(md.mesh.y)); point2=find(md.mesh.x==max(md.mesh.x));
	costheta=(md.mesh.x(point2)-md.mesh.x(point1))/sqrt((md.mesh.x(point2)-md.mesh.x(point1)).^2+(md.mesh.y(point2)-md.mesh.y(point1)).^2);
	sintheta=(md.mesh.y(point2)-md.mesh.y(point1))/sqrt((md.mesh.x(point2)-md.mesh.x(point1)).^2+(md.mesh.y(point2)-md.mesh.y(point1)).^2);
	md.stressbalance.referential(:,1:3)=repmat([costheta,sintheta,0],md.mesh.numberofvertices,1);
	md.stressbalance.referential(:,4:6)=repmat([-sintheta,costheta,0],md.mesh.numberofvertices,1);

	md.basalforcings.groundedice_melting_rate=zeros(md.mesh.numberofvertices,1);
	md.basalforcings.floatingice_melting_rate=zeros(md.mesh.numberofvertices,1);
	md.thermal.spctemperature=[md.initialization.temperature; 1]; %impose observed temperature on surface
	md.basalforcings.geothermalflux=.064*ones(md.mesh.numberofvertices,1);

	%flow equations:
	md=setflowequation(md,'SSA','all');

	savemodel(org,md);
end
% }}}
% {{{ RunUncoupledISSM:
if perform(org,'RunUncoupledISSM'),

	loaddata(org,'Parameters');
	loaddata(org,'ParameterizeIce');

	%timestepping:
	md.timestepping.final_time=100;
	md.timestepping.time_step=0.5;
	md.transient.isgroundingline=0;
	md.transient.isthermal=0;
	md.groundingline.migration='SubelementMigration';
	md.groundingline.melt_interpolation='NoMeltOnPartiallyFloating';
	md.groundingline.friction_interpolation='SubelementFriction2';
	md.masstransport.stabilization=1;

	md.cluster=generic('name',oshostname(),'np',2);
	md=solve(md,'Transient');

	savemodel(org,md);

	plotmodel(md,'data',md.results.TransientSolution(end).Vel,'data',md.results.TransientSolution(end).Thickness)
end
% }}}

%Run MITgcm/ISSM
% {{{ RunCoupledMITgcmISSM:
if perform(org,'RunCoupledMITgcmISSM'),

    loaddata(org,'Parameters');
    loaddata(org,'Bathymetry');
    loaddata(org,'IceSheetGeometry');
    loaddata(org,'ParameterizeIce');
    md=loadmodel(org,'RunUncoupledISSM');

    % {{{ prepare ISSM: start from the steady-state

    md.geometry.base=md.results.TransientSolution(end).Base;
    md.geometry.surface=md.results.TransientSolution(end).Surface;
    md.geometry.thickness=md.results.TransientSolution(end).Thickness;
    md.initialization.vx=md.results.TransientSolution(end).Vx;
    md.initialization.vy=md.results.TransientSolution(end).Vy;
    md.initialization.vel=md.results.TransientSolution(end).Vel;
    md.initialization.pressure=md.results.TransientSolution(end).Pressure;
    md.transient.isoceancoupling=2;
    md.transient.isgroundingline=0;
    md.masstransport.requested_outputs={'default','BasalforcingsFloatingiceMeltingRate'};

    % }}}
    % {{{ prepare MITgcm
    % rename previous run directory and create new one
    if exist ('run.old')
        !\rm -rf run.old
    end
    if exist ('run')
        !\mv run run.old
    end
    !\mkdir run
    !\cp ../MITgcm/build/mitgcmuv run
    !\cp ../MITgcm/input_4003/* run
    
    % {{{ Construct MITgcm binary input files
    namF='run/bathy_flat.bin';
    fid=fopen(namF,'w','b'); fwrite(fid,bathymetry,prec);fclose(fid);
    
    draft=md.results.TransientSolution(end).Base;
    draft=reshape(draft,[Nx+1,Ny+1]);
    pos=ones(Nx+1,Ny+1);
    pos(find(md.mask.ice_levelset>0))=0;
    draft=draft.*pos;
    draft=draft(1:Nx,1:Ny)+draft(2:end,1:Ny)+draft(1:Nx,2:end)+draft(2:end,2:end);
    pos=pos(1:Nx,1:Ny)+pos(2:end,1:Ny)+pos(1:Nx,2:end)+pos(2:end,2:end);
    draft(find(pos))=draft(find(pos))./pos(find(pos));
    namF='run/shelficeTopo.Lin.bin';
    fid=fopen(namF,'w','b'); fwrite(fid,draft,prec);fclose(fid);
    
    del_T = (T_bot - T_sfc)/(59*delZ);
    tref=zeros(1,Nz);
    for k = 1:Nz;
        tref(k) = T_sfc + del_T*((k-20)*delZ);
        tref(k)= max(T_sfc,min(tref(k),T_bot));
    end
    namF='run/temp_obc.bin';
    tref=[tref; tref; tref];
    fid=fopen(namF,'w','b'); fwrite(fid,tref,prec);fclose(fid);

    del_S = (S_bot - S_sfc)/(59*delZ);
    sref=zeros(1,Nz);
    for k = 1:Nz;
        sref(k) = S_sfc + del_S*((k-20)*delZ);
        sref(k)= max(S_sfc,min(sref(k),S_bot));
    end
    namF='run/salt_obc.bin';
    sref=[sref; sref; sref];
    fid=fopen(namF,'w','b'); fwrite(fid,sref,prec);fclose(fid);

    zax=[1:Nz];
    v1=2.5e-2;
    var=1+Nz-2*zax; var=var/(Nz-1);
    vobc=v1*var;
    namF='run/vVel_obc.bin';
    vobc=[vobc; vobc; vobc];
    fid=fopen(namF,'w','b'); fwrite(fid,vobc,prec);fclose(fid);

    var=zeros(Nx,Ny,Nz);
    for i=1:Nx, for j=1:Ny
            var(i,j,:)=tref(1,:);
        end, end
    namF='run/temp_ini.bin';
    fid=fopen(namF,'w','b'); fwrite(fid,var,prec);fclose(fid);

    for i=1:Nx, for j=1:Ny
            var(i,j,:)=sref(1,:);
        end, end
    namF='run/salt_ini.bin';
    fid=fopen(namF,'w','b'); fwrite(fid,var,prec);fclose(fid);
    % }}}
    % }}}

    md.timestepping.coupling_time=1/24/365;          % 1 hour in decimal years
    md.timestepping.time_step=1/24/365;              % 1 hour in decimal years
    md.timestepping.final_time=23/24/365;            % hour 23 in decimal years
    md.cluster.npocean=nPx*nPy;
    md.cluster.np=2;
    md.cluster.executionpath=[pwd '/run'];
    md.transient.requested_outputs={'default','MaskOceanLevelset'};

    md=solveiceocean(md,'Transient','runtimename',false);

%eval(['!mpiexec -np ' int2str(md.cluster.np) ' ' md.cluster.codepath '/issm_ocean.exe TransientSolution ' pwd ' ' md.miscellaneous.name ' ']);
%eval(['!mpiexec -np ' int2str(md.cluster.np) ' ' md.cluster.codepath '/issm_ocean.exe TransientSolution ' pwd ' ' md.miscellaneous.name ' : -np ' int2str(nPx*nPy) ' ./mitgcmuv']);
end
% }}}
