function output = gamv(x,y,data,varargin)
%GAMV - use gslib for Kriging
%
%   Usage:
%      output = gamv(x,y,data,varargin)

options=pairoptions(varargin{:});

nlag = getfieldvalue(options,'nlag', 20);
dlag = getfieldvalue(options,'dlag', 1000);

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
fid=fopen('gamv.par','w');
fprintf(fid,'\t\t\t\t%s\n','Parameters for GAMV');
fprintf(fid,'\t\t\t\t%s\n','*******************');
fprintf(fid,'\n');
fprintf(fid,'%s\n','START OF PARAMETERS:');
fprintf(fid,'%-30s %s\n','./cluster.dat'              ,'\file with data');
fprintf(fid,'%-30s %s\n','1 2 0'                      ,'\columns for X, Y, Z coordinates');
fprintf(fid,'%-30s %s\n','1 3  '                      ,'\number of variables, column number');
fprintf(fid,'%-30s %s\n','-1.0e21 1.0e21'             ,'\trimming limits');
fprintf(fid,'%-30s %s\n','gamv.out'                   ,'\file for variogram output');
fprintf(fid,'%-30s %s\n',num2str(nlag,'%i')           ,'\number of lags');
fprintf(fid,'%-30s %s\n',num2str(dlag,'%g')           ,'\lag separation distance');
fprintf(fid,'%-30s %s\n',num2str(dlag/2,'%g')         ,'\lag tolerance');
fprintf(fid,'%-30s %s\n','3'                          ,'\number of directions');
fprintf(fid,'%-30s %s\n','0.0 90.0 50.0 0.0 90.0 50.0','\azm, atol, bandh, dip, dtol, bandv');
fprintf(fid,'%-30s %s\n','0.0 22.5 25.0 0.0 22.5 25.0','\azm, atol, bandh, dip, dtol, bandv');
fprintf(fid,'%-30s %s\n','90. 22.5 25.0 0.0 22.5 25.0','\azm, atol, bandh, dip, dtol, bandv');
fprintf(fid,'%-30s %s\n','0'                          ,'\standardize sill? (0=no, 1=yes)');
fprintf(fid,'%-30s %s\n','2'                          ,'\number of variograms');
fprintf(fid,'%-30s %s\n','1 1 1'                      ,'\tail var., head vars., variogram type');
fprintf(fid,'%-30s %s\n','1 1 3'                      ,'\tail var., head vars., variogram type');
fclose(fid);

%Call gamv
system([issmdir() '/externalpackages/gslib/install/gamv gamv.par']);
delete('gamv.par');

%Read output
output   = struct('Semivariogram',[],'Covariance',[]);
counter1 = 1;
counter2 = 1;
fid=fopen('gamv.out','r');
while (~feof(fid)),
	A=fscanf(fid,'%s',1);
	if strcmp(A,'Covariance');
		A=fscanf(fid,'%s',4); %Read tail:Data head:Data direction  2
		output(counter1).Covariance=fscanf(fid,'%i %g %g %i %g %g',[6 nlag+2])';
		counter1=counter1+1;
	elseif strcmp(A,'Semivariogram'),
		A=fscanf(fid,'%s',4); %Read tail:Data head:Data direction  2
		output(counter2).Semivariogram=fscanf(fid,'%i %g %g %i %g %g',[6 nlag+2])';
		counter2=counter2+1;
	else
		%do nothing
	end
end
fclose(fid);
delete('gamv.out')
