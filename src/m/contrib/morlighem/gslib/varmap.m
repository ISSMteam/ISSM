function output = varmap(x,y,data,varargin)
%VARMAP - use gslib for Kriging
%
%   Usage:
%      output = varmap(x,y,data,varargin)

options=pairoptions(varargin{:});

nxlag = getfieldvalue(options,'nxlag', 20);
nylag = getfieldvalue(options,'nylag', 20);
dxlag = getfieldvalue(options,'dxlag', 1000);
dylag = getfieldvalue(options,'dylag', 1000);

%Write data file
fid=fopen('cluster.dat','w');
fprintf(fid,'%s\n','Data file');
fprintf(fid,'%i\n',3);
fprintf(fid,'%s\n','Xlocation');
fprintf(fid,'%s\n','Ylocation');
fprintf(fid,'%s\n','Data');
fprintf(fid,'%g %g %g\n',[x y data]');
fclose(fid);

%Write parameter file
fid=fopen('varmap.par','w');
fprintf(fid,'\t\t\t\t%s\n','Parameters for GAMV');
fprintf(fid,'\t\t\t\t%s\n','*******************');
fprintf(fid,'\n');
fprintf(fid,'%s\n','START OF PARAMETERS:');
fprintf(fid,'%-30s %s\n','./cluster.dat'              ,'\file with data');
fprintf(fid,'%-30s %s\n','1 3  '                      ,'\number of variables, column number');
fprintf(fid,'%-30s %s\n','-1.0e21 1.0e21'             ,'\trimming limits');
fprintf(fid,'%-30s %s\n','0    '                      ,'\1=regular grid, 0=scattered values');
fprintf(fid,'%-30s %s\n','50 50 1'                    ,'\if =1: nx, ny, nz');
fprintf(fid,'%-30s %s\n','1.0 1.0 1.0'                ,'\       xsiz, ysiz, zsiz if igrid=1');
fprintf(fid,'%-30s %s\n','1 2 0'                      ,'\if =0: columns for x, y and z coordinates');
fprintf(fid,'%-30s %s\n','varmap.out'                 ,'\file for variogram output');
fprintf(fid,'%-30s %s\n',num2str([nxlag nylag 0],'%i '),'\nxlag, nylag, nzlag');
fprintf(fid,'%-30s %s\n',num2str([dxlag dylag 1],'%g %g %i'),'\dxlag, dylag, dzlag');
fprintf(fid,'%-30s %s\n','5'                          ,'\minimum number of pairs');
fprintf(fid,'%-30s %s\n','0'                          ,'\standardize sill? (0=no, 1=yes)');
fprintf(fid,'%-30s %s\n','1'                          ,'\number of variograms');
fprintf(fid,'%-30s %s\n','1 1 1'                      ,'\tail, head, variogram type');
fclose(fid);

%Call varmap
system([issmdir() '/externalpackages/gslib/install/varmap varmap.par']);
delete('varmap.par');

%Read output
fid=fopen('varmap.out','r');
A = textscan(fid,'%f %f %f %f %f %f','headerlines',8);
fclose(fid);
delete('varmap.out')
output = reshape(A{1},[2*nxlag+1 2*nylag+1]);
