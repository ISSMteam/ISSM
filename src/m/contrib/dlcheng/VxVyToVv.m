name='CWGreenland.'; %new Central West Greenland project to do AD from Store to Upernavik.
org=organizer('repository','./Models','prefix',name,'steps',steps);
NSIDCName='./Data/Greenland';
WestCoastName='./Data/Rink/';
modeldatapath='/Users/larour/ModelData/';
modeldata2path='/Users/larour/ModelData2/';
MARName=[modeldata2path 'MARv3.11/ERA/CW_Greenland/'];
qgispath='/Users/larour/issm-jpl/proj-group/qgis/';
greenlandshppath=[qgispath '/GreenlandStateEstimate/'];
GIMP='/Users/larour/ModelData/HowatDEMGreenland2012/HowatDEMGreenland2012_90m_smooth.mat';

%Path to shapefiles: 
shppath=[qgispath '/Greenland/30km/'];
shppathslr=[qgispath '/Slr/'];

if perform(org,'GenerateLandsatVV'), % {{{
    rootname = [modeldatapath '/VelHowat'];
    infos = dir([rootname '/*/*/*_vx_*.tif']);
    for i=1:length(infos),
        filename = infos(i).name; %Example: OPT_W70.55N_1986-05_vx_v02.1.tif
        path = [infos(i).folder '/' infos(i).name];

        vx_path = path;
        vy_path = strrep(vx_path, '_vx_', '_vy_');
        vv_path = strrep(vx_path, '_vx_', '_vv_');

        [vx,Rx] = readgeoraster(vx_path);
        [vy,Ry] = readgeoraster(vy_path);
        vx(vx==-99999)=nan;
        vy(vy==-99999)=nan;

        vel = sqrt(vx.^2+vy.^2);
        disp(vv_path);
        geotiffwrite(vv_path, vel, Rx, 'CoordRefSysCode', 3413);
    end
end % }}}

