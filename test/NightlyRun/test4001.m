%Test Name: IceOcean
%ISSM/MITgcm coupled set-up
%
%Script control parameters
steps=1:11;
final_time=1;

!rm -rf ${ISSM_DIR}/test/MITgcm/install
!rm -rf ${ISSM_DIR}/test/MITgcm/build/*
!rm -f ${ISSM_DIR}/test/MITgcm/code/SIZE.h
!rm -rf Models

%Organizer
!mkdir Models
org=organizer('repository','Models','prefix','IceOcean.','steps',steps);

presentdirectory=pwd;

% {{{ Parameters:
if perform(org,'Parameters'),
    Nx=20; %number of longitude cells
    Ny=40; %number of latitude cells
    Nz=30; %number of MITgcm vertical cells
    nPx=2; %number of MITgcm processes to use in x direction
    nPy=4; %number of MITgcm processes to use in y direction
    xgOrigin=0; %origin of longitude
    ygOrigin=-80; %origin of latitude
    dLong=.25; %longitude grid spacing
    dLat=.05; %latitude grid spacing
    delZ=30; %thickness of vertical levels
    icefront_position_ratio=.75;
    ice_thickness=100;
    rho_ice=917;
    rho_water=1028.14;
    di=rho_ice/rho_water;

    % MITgcm initial and lateral boundary conditions
    iniSalt  = 34.4; % initial salinity (PSU)
    iniTheta = -1.9; % initial potential temperature (deg C)
    obcSalt  = 34.4; % open boundary salinity (PSU)
    obcTheta =  1.0; % open boundary potential temperature (deg C)
    mlDepth  = 120.; % mixed layer depth (m)
    mlSalt   = 33.4; % open boundary salinity (PSU)
    mlTheta  = -1.9; % open boundary potential temperature (deg C)
    obcUvel  = -0.1; % open boundary velocity (m/s)

    MITgcmDeltaT=600; % MITgcm time step in seconds
    y2s=31536000; % year to seconds conversion, i.e., seconds per year

    % start_time, final_time, and time_step
    start_time=0; % in decimal years
    time_step=1/12; % coupling interval in decimal years
    async_step_MITgcm_multiplier=1/30; % used to reduce run time for MItgcm

    % bedrock/bathymetry
    hmax=1000;
    trough_depth=200;
    deltah=300;
    sea_level=1095;

    % issm settings:
    numlayers=10;

    savedata(org, Nx, Ny, nPx, nPy, Nz, dLong, dLat, delZ, xgOrigin, ...
        ygOrigin, icefront_position_ratio, ice_thickness, rho_ice, ...
        rho_water, di, hmax, trough_depth, deltah, sea_level, ...
        iniSalt, iniTheta, obcSalt, obcTheta, mlDepth, mlSalt, ...
        mlTheta, obcUvel, start_time, time_step, MITgcmDeltaT, y2s,...
        numlayers,async_step_MITgcm_multiplier);
end
% }}}
% {{{ Bathymetry:
if perform(org,'Bathymetry'),

    loaddata(org,'Parameters');
    %create lat,long
    lat=(ygOrigin+dLat/2):dLat:(ygOrigin+Ny*dLat);
    long=(xgOrigin+dLong/2):dLong:(xgOrigin+Nx*dLong);
    [lat long]=meshgrid(lat,long);

    longmin=min(long(:));
    longmax=max(long(:));
    latmin=min(lat(:));
    latmax=max(lat(:));

    %create bedrock/bathymetry:
    bedrock=zeros(Nx,Ny);
    bedrock=hmax-deltah*tanh(pi*(2*(lat-latmin)./(latmax-latmin)-1))+ ...
            trough_depth*cos(2*pi*long./(longmax-longmin));

    %save bathymetry file for MITgcm
    bathymetry=bedrock-sea_level;
    savedata(org,lat,long,bathymetry);

end
% }}}
% {{{ IceSheetGeometry:
if perform(org,'IceSheetGeometry'),

    loaddata(org,'Parameters');
    loaddata(org,'Bathymetry');
    latmin=min(lat(:));
    latmax=max(lat(:));

    %put ice_thickness constant layer of ice over the bathymetry, unless it floats:
    s=size(bathymetry);
    thickness=ice_thickness*ones(s);

    %figure out ice shelf:
    pos=find(-di*thickness>bathymetry);
    iceshelf_mask=zeros(s);
    iceshelf_mask(pos)=1;

    ice_mask=ones(s);
    pos=find((lat-latmin)/(latmax-latmin)>(icefront_position_ratio));
    ice_mask(pos)=0;
    iceshelf_mask(pos)=0;

    %compute draft of ice shelf:
    draft=bathymetry;
    pos=find(iceshelf_mask);
    draft(pos)=-di*thickness(pos);
    pos=find(~ice_mask);
    draft(pos)=0;

    savedata(org,ice_mask,iceshelf_mask,draft,thickness);
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

    %load data:
    loaddata(org,'Parameters');

    %specify computational grid in SIZE.h
    if ~exist('../MITgcm/code/SIZE.h')
        fidi=fopen('../MITgcm/code/SIZE.h.bak','r');
        fido=fopen('../MITgcm/code/SIZE.h','w');
        tline = fgetl(fidi);
        fprintf(fido,'%s\n',tline);
        while 1
            tline = fgetl(fidi);
            if ~ischar(tline), break, end
            %do the change here:
            if strcmpi(tline,'     &           sNx =  20,'),
                fprintf(fido,'%s%i%s\n','     &           sNx =  ',round(Nx/nPx),',');
                continue;
            end
            if strcmpi(tline,'     &           sNy =  20,'),
                fprintf(fido,'%s%i%s\n','     &           sNy =  ',round(Ny/nPy),',');
                continue;
            end
            if strcmpi(tline,'     &           nPx =   1,'),
                fprintf(fido,'%s%i%s\n','     &           nPx = ',nPx,',');
                continue;
            end
            if strcmpi(tline,'     &           nPy =   2,'),
                fprintf(fido,'%s%i%s\n','     &           nPy = ',nPy,',');
                continue;
            end
            fprintf(fido,'%s\n',tline);
        end
        %close  files
        fclose(fidi);
        fclose(fido);
    end

    system(['../MITgcm/build.sh generic ' pwd '/../MITgcm']);
end
% }}}
addpath(recursivepath([pwd '/../MITgcm']));
% {{{ RunUncoupledMITgcm:
if perform(org,'RunUncoupledMITgcm'),

    %load data:
    loaddata(org,'Parameters');
    loaddata(org,'Bathymetry');
    loaddata(org,'IceSheetGeometry');
     endtime = round(MITgcmDeltaT * ...
         floor(time_step*y2s*async_step_MITgcm_multiplier/MITgcmDeltaT));

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
    !\cp ../MITgcm/input/* run
    !\cp ../MITgcm/input/eedata_uncoupled run/eedata

    %load data:
    loaddata(org,'Parameters');

    % initial salinity
    S=iniSalt*ones(Nx,Ny,Nz);
    writebin('run/Salt.bin',S);

    % initial temperature
    T=iniTheta*ones(Nx,Ny,Nz);
    writebin('run/Theta.bin',T);

    % initial velocity
    Z=zeros(Nx,Ny,Nz);
    writebin('run/Uvel.bin',Z);
    writebin('run/Vvel.bin',Z);

    % initial sea surface height
    Z=zeros(Nx,Ny);
    writebin('run/Etan.bin',Z);

    % salinity boundary conditions
    S=obcSalt*ones(Ny,Nz);
    thk=delZ*ones(Nz,1);
    bot=cumsum(thk);
    ik=find(bot<=mlDepth);
    S(:,ik)=mlSalt;
    writebin('run/OBs.bin',S);

    % temperature boundary conditions
    T=obcTheta*ones(Ny,Nz);
    T(:,ik)=mlTheta;
    writebin('run/OBt.bin',T);

    % zonal velocity boundary conditions
    U=obcUvel*ones(Ny,Nz);
    writebin('run/OBu.bin',U);

    % zero boundary conditions
    Z=zeros(Ny,Nz);
    writebin('run/zeros.bin',Z);

    % build parameter file data.obcs
    fidi=fopen('../MITgcm/input/data.obcs','r');
    fido=fopen('run/data.obcs','w');
    tline = fgetl(fidi);
    fprintf(fido,'%s\n',tline);
    while 1
        tline = fgetl(fidi);
        if ~ischar(tline), break, end
        %do the change here:
        if strcmpi(tline,' OB_Iwest = 40*1,'),
            fprintf(fido,'%s%i%s\n',' OB_Iwest = ',Ny,'*1,');
            continue;
        end
        if strcmpi(tline,' OB_Ieast = 40*-1,'),
            fprintf(fido,'%s%i%s\n',' OB_Ieast = ',Ny,'*-1,');
            continue;
        end
        fprintf(fido,'%s\n',tline);
    end
    %close  files
    fclose(fidi);
    fclose(fido);

    %save bathymetry and bedrock in run directory
    writebin('run/bathymetry.bin',bathymetry);
    writebin('run/icetopo.bin',draft);
    % }}}

    %start looping:
    for t=start_time:time_step:final_time,
        disp(['Year: ' num2str(t)])
        % {{{ generate MITgcm parameter file data
        fidi=fopen('../MITgcm/input/data','r');
        fido=fopen('run/data','w');
        tline = fgetl(fidi);
        fprintf(fido,'%s\n',tline);
        while 1
            tline = fgetl(fidi);
            if ~ischar(tline), break, end
            %do the change here:
            if strcmpi(tline,' xgOrigin = 0.0,'),
                fprintf(fido,'%s%i%s\n',' xgOrigin = ',xgOrigin,',');
                continue;
            end
            if strcmpi(tline,' ygOrigin = -80.0,'),
                fprintf(fido,'%s%i%s\n',' ygOrigin = ',ygOrigin,',');
                continue;
            end
            if strcmpi(tline,' delX = 20*0.25,'),
                fprintf(fido,'%s%i*%g%s\n',' delX = ',Nx,dLong,',');
                continue;
            end
            if strcmpi(tline,' delY = 20*0.25,'),
                fprintf(fido,'%s%i*%g%s\n',' delY = ',Ny,dLat,',');
                continue;
            end
            if strcmpi(tline,' delZ = 30*30.0,'),
                fprintf(fido,'%s%i*%g%s\n',' delZ = ',Nz,delZ,',');
                continue;
            end
            if strcmpi(tline,' endTime=2592000.,'),
                fprintf(fido,'%s%i%s\n',' endTime= ',endtime,',');
                continue;
            end
            if strcmpi(tline,' deltaT=1200.0,'),
                fprintf(fido,'%s%i%s\n',' deltaT= ',MITgcmDeltaT,',');
                continue;
            end
            if strcmpi(tline,' pChkptFreq=2592000.,'),
                fprintf(fido,'%s%i%s\n',' pChkptFreq= ',endtime,',');
                continue;
            end
            if strcmpi(tline,' taveFreq=2592000.,'),
                fprintf(fido,'%s%i%s\n',' taveFreq= ',endtime,',');
                continue;
            end
            if strcmpi(tline,' rhoConst=1030.,'),
                fprintf(fido,'%s%i%s\n',' rhoConst= ',rho_water,',');
                continue;
            end
            if strcmpi(tline,' rhoNil=1030.,'),
                fprintf(fido,'%s%i%s\n',' rhoNil= ',rho_water,',');
                continue;
            end
            fprintf(fido,'%s\n',tline);
        end
        %close  files
        fclose(fidi);
        fclose(fido);
        % }}}
        % {{{ generate initial MITgcm conditions

        ds=round(endtime/MITgcmDeltaT);
        if t>start_time
            % Read pickup file
            fnm=['run/pickup.' myint2str(ds,10) '.data'];
            U=readbin(fnm,[Nx Ny Nz],1,'real*8',0);
            V=readbin(fnm,[Nx Ny Nz],1,'real*8',1);
            T=readbin(fnm,[Nx Ny Nz],1,'real*8',2);
            S=readbin(fnm,[Nx Ny Nz],1,'real*8',3);
            E=readbin(fnm,[Nx Ny],1,'real*8',8*Nz);
            writebin('run/Salt.bin' ,S);
            writebin('run/Theta.bin',T);
            writebin('run/Uvel.bin' ,U);
            writebin('run/Vvel.bin' ,V);
            writebin('run/Etan.bin' ,E);
        end

        % }}}
        % {{{ system call to run MITgcm
        cd run
        eval(['!mpirun -np ' int2str(nPx*nPy) ' ./mitgcmuv']);
        ts=round((t+time_step)*y2s/MITgcmDeltaT);
        eval(['!\mv STDERR.0000 STDERR_' myint2str(ts,10) '.data'])
        eval(['!\mv STDOUT.0000 STDOUT_' myint2str(ts,10) '.data'])
        eval(['!\cp hFacC.data hFacC_' myint2str(ts,10) '.data'])
        eval(['!\cp icetopo.bin icetopo_' myint2str(ts,10) '.data'])
        for fld={'S','T','U','V','Eta', ...
                 'SHICE_heatFluxtave','SHICE_fwFluxtave'}
            eval(['!\mv ' fld{1} '.' myint2str(ds,10) '.data ' ...
                  fld{1} '_' myint2str(ts,10) '.data'])
        end
        cd ..
        % }}}
    end
end
% }}}

%Configure ISSM
% {{{ CreateMesh:
if perform(org,'CreateMesh'),

    loaddata(org,'Parameters');
    loaddata(org,'Bathymetry');
    loaddata(org,'IceSheetGeometry');

    %create model:
    md=model();

    %Grab lat,long from MITgcm:
    lat=lat(:);
    long=long(:);

    %project lat,long:
    [x,y]=ll2xy(lat,long,-1);

    index=[];
    %  C  D
    %  A  B
    for j=1:Ny-1,
        for i=1:Nx-1,
            A=(j-1)*Nx+i;
            B=(j-1)*Nx+i+1;
            C=j*Nx+i;
            D=j*Nx+i+1;
            index(end+1,:)=[A B C];
            index(end+1,:)=[C B D];
        end
    end

    %fill mesh and model:
    md=meshconvert(md,index,x,y);
    md.mesh.lat=lat;
    md.mesh.long=long;

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
    bathymetry=bathymetry(:);
    iceshelf_mask=iceshelf_mask(:);
    ice_mask=ice_mask(:);
    thickness=thickness(:);
    draft=draft(:);

    %start filling some of the fields
    md.geometry.bed=bathymetry;
    md.geometry.thickness=thickness;
    md.geometry.base=md.geometry.bed;
    pos=find(iceshelf_mask); md.geometry.base(pos)=draft(pos);
    md.geometry.surface=md.geometry.base+md.geometry.thickness;

    %nothing passes icefront:
    pos=find(~ice_mask);
    md.geometry.thickness(pos)=1;
    md.geometry.surface(pos)=(1-di)*md.geometry.thickness(pos);
    md.geometry.base(pos)=-di*md.geometry.thickness(pos);

    %level sets:
    md.mask.ocean_levelset=-ones(md.mesh.numberofvertices,1);
    md.mask.ice_levelset=ones(md.mesh.numberofvertices,1);

    pos=find(ice_mask); md.mask.ice_levelset(pos)=-1;
    pos=find(~iceshelf_mask & ice_mask); md.mask.ocean_levelset(pos)=1;

    %identify edges of grounded ice:
    ocean_levelset=md.mask.ocean_levelset;
    for i=1:md.mesh.numberofelements,
        m=ocean_levelset(md.mesh.elements(i,:));
        if abs(sum(m))~=3,
            pos=find(m==1); md.mask.ocean_levelset(md.mesh.elements(i,pos))=0;
        end
    end

    %identify edges of ice:
    ice_levelset=md.mask.ice_levelset;
    for i=1:md.mesh.numberofelements,
        m=ice_levelset(md.mesh.elements(i,:));
        if abs(sum(m))~=3,
            pos=find(m==-1); md.mask.ice_levelset(md.mesh.elements(i,pos))=0;
        end
    end

    savemodel(org,md);
end
% }}}
% {{{ ParameterizeIce:
if perform(org,'ParameterizeIce'),

    loaddata(org,'Parameters');
    loaddata(org,'CreateMesh');
    loaddata(org,'MeshGeometry');

    %miscellaneous
    md.miscellaneous.name='test4001';

    %initial velocity:
    md.initialization.vx=zeros(md.mesh.numberofvertices,1);
    md.initialization.vy=zeros(md.mesh.numberofvertices,1);
    md.initialization.vz=zeros(md.mesh.numberofvertices,1);

    %friction:
    md.friction.coefficient=30*ones(md.mesh.numberofvertices,1);
    pos=find(md.mask.ocean_levelset<=0);
    md.friction.coefficient(pos)=0;
    md.friction.p=ones(md.mesh.numberofelements,1);
    md.friction.q=ones(md.mesh.numberofelements,1);

    %temperatures and surface mass balance:
    md.initialization.temperature=(273.15-20)*ones(md.mesh.numberofvertices,1);
    md.initialization.pressure=md.materials.rho_ice*md.constants.g*(md.geometry.surface-md.geometry.base);
    md.smb.mass_balance = [1*ones(md.mesh.numberofvertices,1); 1];

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

    %deal with water:
    pos=find(md.mask.ice_levelset>0);
    md.stressbalance.spcvx(pos)=0;
    md.stressbalance.spcvy(pos)=0;
    md.stressbalance.spcvz(pos)=0;
    md.masstransport.spcthickness(pos)=0;

    %get some flux at the ice divide:
    pos=find(md.mesh.lat==min(md.mesh.lat));
    md.stressbalance.spcvy(pos)=200;

    %deal with boundaries, excluding icefront:
    vertex_on_boundary=zeros(md.mesh.numberofvertices,1);
    vertex_on_boundary(md.mesh.segments(:,1:2))=1;
    pos=find(vertex_on_boundary & md.mask.ocean_levelset<=0);
    md.stressbalance.spcvx(pos)=md.initialization.vx(pos);
    md.stressbalance.spcvy(pos)=md.initialization.vy(pos);
    md.stressbalance.spcvz(pos)=md.initialization.vz(pos);
    md.masstransport.spcthickness(pos)=md.geometry.thickness(pos);

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
    md.timestepping.final_time=final_time;
    md.timestepping.time_step=time_step;
    md.transient.isgroundingline=1;
    md.transient.isthermal=0;
    md.groundingline.migration='SubelementMigration';
    md.groundingline.melt_interpolation='SubelementMelt2';
    md.groundingline.friction_interpolation='SubelementFriction2';

    md.cluster=generic('name',oshostname(),'np',2);
    md=solve(md,'Transient');

    savemodel(org,md);
end
% }}}

%Run MITgcm/ISSM
% {{{ RunCoupledMITgcmISSM:
if perform(org,'RunCoupledMITgcmISSM'),

    %load data:
    loaddata(org,'Parameters');
    loaddata(org,'ParameterizeIce');
    loaddata(org,'Bathymetry');
    loaddata(org,'IceSheetGeometry');
        endtime = round(MITgcmDeltaT * ...
         floor(time_step*y2s*async_step_MITgcm_multiplier/MITgcmDeltaT));

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
        !\cp ../MITgcm/input/* run
        !\cp ../MITgcm/input/eedata_uncoupled run/eedata

        %load data:
        loaddata(org,'Parameters');

        % initial salinity
        S=iniSalt*ones(Nx,Ny,Nz);
        writebin('run/Salt.bin',S);

        % initial temperature
        T=iniTheta*ones(Nx,Ny,Nz);
        writebin('run/Theta.bin',T);

        % initial velocity
        Z=zeros(Nx,Ny,Nz);
        writebin('run/Uvel.bin',Z);
        writebin('run/Vvel.bin',Z);

        % initial sea surface height
        Z=zeros(Nx,Ny);
        writebin('run/Etan.bin',Z);

        % salinity boundary conditions
        S=obcSalt*ones(Ny,Nz);
        thk=delZ*ones(Nz,1);
        bot=cumsum(thk);
        ik=find(bot<=mlDepth);
        S(:,ik)=mlSalt;
        writebin('run/OBs.bin',S);

        % temperature boundary conditions
        T=obcTheta*ones(Ny,Nz);
        T(:,ik)=mlTheta;
        writebin('run/OBt.bin',T);

        % zonal velocity boundary conditions
        U=obcUvel*ones(Ny,Nz);
        writebin('run/OBu.bin',U);

        % zero boundary conditions
        Z=zeros(Ny,Nz);
        writebin('run/zeros.bin',Z);

        % build parameter file data.obcs
        fidi=fopen('../MITgcm/input/data.obcs','r');
        fido=fopen('run/data.obcs','w');
        tline = fgetl(fidi);
        fprintf(fido,'%s\n',tline);
        while 1
            tline = fgetl(fidi);
            if ~ischar(tline), break, end
            %do the change here:
            if strcmpi(tline,' OB_Iwest = 40*1,'),
                fprintf(fido,'%s%i%s\n',' OB_Iwest = ',Ny,'*1,');
                continue;
            end
            if strcmpi(tline,' OB_Ieast = 40*-1,'),
                fprintf(fido,'%s%i%s\n',' OB_Ieast = ',Ny,'*-1,');
                continue;
            end
            fprintf(fido,'%s\n',tline);
        end
        %close  files
        fclose(fidi);
        fclose(fido);

        %save bathymetry in MITgcm run directory
        writebin('run/bathymetry.bin',bathymetry);
        % }}}

    % {{{ ISSM settings:

    setenv('DYLD_LIBRARY_PATH', '/usr/local/gfortran/lib')
    %timestepping:
    md.timestepping.start_time=start_time;
    md.timestepping.final_time=final_time;
    md.timestepping.time_step=time_step;
    md.cluster=generic('name',oshostname(),'np',2);
    md.results.TransientSolution.Base=md.geometry.base;
    md.transient.isgroundingline=1;
    md.transient.isthermal=0;
    md.groundingline.migration='SubelementMigration';
    md.groundingline.melt_interpolation='SubelementMelt2';
    md.groundingline.friction_interpolation='SubelementFriction2';

    % }}}

    %start looping:
    results=md.results;

    for t=start_time:time_step:final_time
            disp(['Year: ' num2str(t)])

        %send draft from ISSM to MITgcm:
        draft=md.results.TransientSolution(end).Base;
        pos=find(md.mask.ice_levelset>0); draft(pos)=0;
            if t>start_time
                old_draft=readbin('run/icetopo.bin',[Nx,Ny]);
            end
            writebin('run/icetopo.bin',draft);

        % {{{ generate MITgcm parameter file data
        fidi=fopen('../MITgcm/input/data','r');
        fido=fopen('run/data','w');
        tline = fgetl(fidi);
        fprintf(fido,'%s\n',tline);
            while 1
                tline = fgetl(fidi);
                if ~ischar(tline), break, end
        %do the change here:
        if strcmpi(tline,' xgOrigin = 0.0,'),
            fprintf(fido,'%s%i%s\n',' xgOrigin = ',xgOrigin,',');
            continue;
        end
        if strcmpi(tline,' ygOrigin = -80.0,'),
            fprintf(fido,'%s%i%s\n',' ygOrigin = ',ygOrigin,',');
            continue;
        end
        if strcmpi(tline,' delX = 20*0.25,'),
            fprintf(fido,'%s%i*%g%s\n',' delX = ',Nx,dLong,',');
            continue;
        end
        if strcmpi(tline,' delY = 20*0.25,'),
            fprintf(fido,'%s%i*%g%s\n',' delY = ',Ny,dLat,',');
            continue;
        end
                if strcmpi(tline,' delZ = 30*30.0,'),
                    fprintf(fido,'%s%i*%g%s\n',' delZ = ',Nz,delZ,',');
                    continue;
                end
        if strcmpi(tline,' endTime=2592000.,'),
            fprintf(fido,'%s%i%s\n',' endTime= ',endtime,',');
            continue;
        end
        if strcmpi(tline,' deltaT=1200.0,'),
            fprintf(fido,'%s%i%s\n',' deltaT= ',MITgcmDeltaT,',');
            continue;
        end
        if strcmpi(tline,' pChkptFreq=2592000.,'),
            fprintf(fido,'%s%i%s\n',' pChkptFreq= ',endtime,',');
            continue;
        end
        if strcmpi(tline,' taveFreq=2592000.,'),
            fprintf(fido,'%s%i%s\n',' taveFreq= ',endtime,',');
            continue;
        end
                if strcmpi(tline,' rhoConst=1030.,'),
                    fprintf(fido,'%s%i%s\n',' rhoConst= ',rho_water,',');
                    continue;
                end
                if strcmpi(tline,' rhoNil=1030.,'),
                    fprintf(fido,'%s%i%s\n',' rhoNil= ',rho_water,',');
                    continue;
                end
        fprintf(fido,'%s\n',tline);
        end
        %close  files
        fclose(fidi);
        fclose(fido);
        % }}}

        % {{{ generate initial MITgcm conditions
            ds=round(endtime/MITgcmDeltaT);
            if t>start_time
                % Read pickup file
                fnm=['run/pickup.' myint2str(ds,10) '.data'];
                U=readbin(fnm,[Nx Ny Nz],1,'real*8',0);
                V=readbin(fnm,[Nx Ny Nz],1,'real*8',1);
                T=readbin(fnm,[Nx Ny Nz],1,'real*8',2);
                S=readbin(fnm,[Nx Ny Nz],1,'real*8',3);
                E=readbin(fnm,[Nx Ny],1,'real*8',8*Nz);

                % find indices of locations where ice shelf retreated
                h=readbin('run/hFacC.data',[Nx Ny Nz]);
                msk=sum(h,3);
                msk(find(msk))=1;
                [iw jw]=find(msk); % horizontal indices where there is water
                tmp=reshape(draft,[Nx,Ny])-old_draft;
                tmp(find(tmp<0))=0;
                [im jm]=find(tmp); % horizontal indices where there is melt

                % Extrapolate T/S to locations where ice shelf retreated
                for i=1:length(im)

                    % first try vertical extrapolation
                    in=find(h(im(i),jm(i),:));
                    if length(in)>0;
                        S(im(i),jm(i),1:min(in)  ) = S(im(i),jm(i),min(in));
                        T(im(i),jm(i),1:min(in)  ) = T(im(i),jm(i),min(in));
                        continue
                    end

                    % if not succesful, use closest neighbor horizontal extrapolation
                    [y c]=min((iw-im(i)).^2+(jw-jm(i)).^2);
                    salt=squeeze(S(iw(c),jw(c),:)); % salinity profile of closest neighbor
                    temp=squeeze(T(iw(c),jw(c),:)); % salinity profile of closest neighbor
                    in=find(h(iw(c),jw(c),:));
                    salt(1:min(in))=salt(min(in));
                    temp(1:min(in))=temp(min(in));
                    salt(max(in):end)=salt(max(in));
                    temp(max(in):end)=temp(max(in));
                    S(im(i),jm(i),:)=salt;
                    T(im(i),jm(i),:)=temp;
                end

                % Write initial conditions
                writebin('run/Salt.bin' ,S);
                writebin('run/Theta.bin',T);
                writebin('run/Uvel.bin' ,U);
                writebin('run/Vvel.bin' ,V);
                writebin('run/Etan.bin' ,E);
            end
            % }}}

            % {{{ system call to run MITgcm
            cd run
            eval(['!mpirun -np ' int2str(nPx*nPy) ' ./mitgcmuv']);
            ts=round((t+time_step)*y2s/MITgcmDeltaT);
            eval(['!\mv STDERR.0000 STDERR_' myint2str(ts,10) '.data'])
            eval(['!\mv STDOUT.0000 STDOUT_' myint2str(ts,10) '.data'])
            eval(['!\cp hFacC.data hFacC_' myint2str(ts,10) '.data'])
            eval(['!\cp icetopo.bin icetopo_' myint2str(ts,10) '.data'])
            for fld={'S','T','U','V','Eta', ...
                     'SHICE_heatFluxtave','SHICE_fwFluxtave'}
                eval(['!\mv ' fld{1} '.' myint2str(ds,10) '.data ' ...
                      fld{1} '_' myint2str(ts,10) '.data'])
            end
            cd ..
            % }}}

        %get melting rates from MITgcm
        %upward fresh water flux (kg/m^2/s):
        fnm=['run/SHICE_fwFluxtave_' myint2str(ts,10) '.data'];
        melting_rate=readbin(fnm,[Nx Ny]);

        %send averaged melting rate to ISSM
        %downward fresh water flux (m/y):
        melting_rate=-melting_rate(:)*y2s/rho_ice;
         md.basalforcings.floatingice_melting_rate=melting_rate;

        % {{{ run ISSM and recover results

        md.timestepping.start_time=t;
        md.timestepping.final_time=t+time_step;;
        md=solve(md,'Transient');

        base=md.results.TransientSolution(end).Base;
        thickness=md.results.TransientSolution(end).Thickness;
        md.geometry.base=base;
        md.geometry.thickness=thickness;
        md.geometry.surface=md.geometry.base+md.geometry.thickness;
        md.initialization.vx=md.results.TransientSolution(end).Vx;
        md.initialization.vy=md.results.TransientSolution(end).Vy;
        md.initialization.vel=md.results.TransientSolution(end).Vel;
        md.initialization.pressure=md.results.TransientSolution(end).Pressure;
        md.mask.ocean_levelset=md.results.TransientSolution(end).MaskOceanLevelset;
        md.results.TransientSolution(end).FloatingiceMeltingRate=md.basalforcings.floatingice_melting_rate;

        %save these results in the model, otherwise, they'll be wiped out
        results(end+1)=md.results;

        % }}}

    end

    md.results=results;
    savemodel(org,md);
end
% }}}

%Fields and tolerances to track changes
fnm=['run/SHICE_fwFluxtave_0000004380.data'];
melting_rate_1=readbin(fnm,[Nx Ny]);
fnm=['run/SHICE_fwFluxtave_0000008760.data'];
melting_rate_2=readbin(fnm,[Nx Ny]);
fnm=['run/SHICE_fwFluxtave_0000013140.data'];
melting_rate_3=readbin(fnm,[Nx Ny]);
fnm=['run/SHICE_fwFluxtave_0000017520.data'];
melting_rate_4=readbin(fnm,[Nx Ny]);
field_names     ={'Base1','Melting1','Vx2','Vy2','Thickness2','Base2','MaskOceanLevelset2','FloatingiceMeltingRate2',...
    'Melting2','Vx3','Vy3','Thickness3','Base3','MaskOceanLevelset3','FloatingiceMeltingRate3',...
    'Melting3','Vx4','Vy4','Thickness4','Base4','MaskOceanLevelset4','FloatingiceMeltingRate4','Melting4'};
field_tolerances={2e-13,1e-13,...
    8e-06,7e-06,2e-07,2e-08,3e-08,1e-13,1e-13,...
    8e-06,7e-06,4e-07,3e-08,5e-08,1e-13,6e-08,...
    8e-06,7e-06,5e-07,4e-08,8e-08,6e-08,3e-07};
field_values={...
    (md.results(1).TransientSolution(end).Base),...
    (melting_rate_1(:)),...
    (md.results(2).TransientSolution(end).Vx),...
    (md.results(2).TransientSolution(end).Vy),...
    (md.results(2).TransientSolution(end).Thickness),...
    (md.results(2).TransientSolution(end).Base),...
    (md.results(2).TransientSolution(end).MaskOceanLevelset),...
    (md.results(2).TransientSolution(end).FloatingiceMeltingRate),...
    (melting_rate_2(:)),...
    (md.results(3).TransientSolution(end).Vx),...
    (md.results(3).TransientSolution(end).Vy),...
    (md.results(3).TransientSolution(end).Thickness),...
    (md.results(3).TransientSolution(end).Base),...
    (md.results(3).TransientSolution(end).MaskOceanLevelset),...
    (md.results(3).TransientSolution(end).FloatingiceMeltingRate),...
    (melting_rate_3(:)),...
    (md.results(4).TransientSolution(end).Vx),...
    (md.results(4).TransientSolution(end).Vy),...
    (md.results(4).TransientSolution(end).Thickness),...
    (md.results(4).TransientSolution(end).Base),...
    (md.results(4).TransientSolution(end).MaskOceanLevelset),...
    (md.results(4).TransientSolution(end).FloatingiceMeltingRate),...
    (melting_rate_4(:)),...
    };
