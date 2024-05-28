function [B E] = gslib(x,y,data,xmin,ymin,nx,ny,deltax,deltay,varargin)
%GSLIB - use gslib for Kriging
%
%   Usage:
%      output = gslib(x,y,data,varargin)

%process options
options = pairoptions(varargin{:});

%Variogram
nugget= getfieldvalue(options,'nugget',10);
sill  = getfieldvalue(options,'sill',164);
range = getfieldvalue(options,'range',25763);

%Kriging options
mindata = getfieldvalue(options,'mindata',1);
maxdata = getfieldvalue(options,'maxdata',50);
maxsearchradius = getfieldvalue(options,'searchrange',50000);

%Some intermediaries (Convert to gslib's parameters);
c = (sill-nugget);
a = sqrt(3)*range;

%Write data file
fid=fopen('cluster.dat','w');
fprintf(fid,'%s\n','Data file');
fprintf(fid,'%i\n',3);
fprintf(fid,'%s\n','Xlocation');
fprintf(fid,'%s\n','Ylocation');
fprintf(fid,'%s\n','Data');
fprintf(fid,'%g %g %g\n',[x y data]');
fclose(fid);

if 0, %GAMV
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
	fprintf(fid,'%-30s %s\n','20'                         ,'\number of lags');
	fprintf(fid,'%-30s %s\n','5.0'                        ,'\lag separation distance');
	fprintf(fid,'%-30s %s\n','3.0'                        ,'\lag tolerance');
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

else, %Kriging KB2D
	%Write parameter file
	fid=fopen('kb2d.par','w');
	fprintf(fid,'\t\t\t\t%s\n','Parameters for KB2D');
	fprintf(fid,'\t\t\t\t%s\n','*******************');
	fprintf(fid,'\n');
	fprintf(fid,'%s\n','START OF PARAMETERS:');
	fprintf(fid,'%-30s %s\n','./cluster.dat'                  ,'\file with data');
	fprintf(fid,'%-30s %s\n','1 2 3'                          ,'\columns for X, Y and variable');
	fprintf(fid,'%-30s %s\n','-1.0e21 1.0e21'                 ,'\trimming limits');
	fprintf(fid,'%-30s %s\n','0'                              ,'\debugging level: 0,1,2,3');
	fprintf(fid,'%-30s %s\n','kb2d.dbg'                       ,'\file for debuggging output');
	fprintf(fid,'%-30s %s\n','kb2d.out'                       ,'\file for kriged output');
	fprintf(fid,'%-30s %s\n',num2str([nx xmin deltax],'%i %10g %6g')  ,'\nx, xmn, xsiz');
	fprintf(fid,'%-30s %s\n',num2str([ny ymin deltay],'%i %10g %6g')  ,'\nx, xmn, xsiz');
	fprintf(fid,'%-30s %s\n','1 1'                            ,'\x and y block discretization');
	fprintf(fid,'%-30s %s\n',num2str([mindata maxdata],'%6g') ,'\min and max data for kriging');
	fprintf(fid,'%-30s %s\n',num2str(maxsearchradius,'%6g')   ,'\max search radius');
	fprintf(fid,'%-30s %s\n','1 2.302'                        ,'\0=SK, 1=OK, (mean if SK)');
	fprintf(fid,'%-30s %s\n',['1 ' num2str(nugget)]           ,'\nst, nugget effect');
	fprintf(fid,'%-30s %s\n',['3 ' num2str([c 0.0 a a],'%10g')],'\it, c, azm, a_max, a_min');
	fclose(fid);

	tic;system([issmdir() '/externalpackages/gslib/install/kb2d kb2d.par']);toc;
	delete('kb2d.par');

	%Read output
	fid=fopen('kb2d.out','r');
	while (~feof(fid)),
		A=fscanf(fid,'%s',1);
		if strcmp(A,'KB2D');
			A=fscanf(fid,'%s',1); %Read output
			params=fscanf(fid,'%i %i %i %i %g %g %g %g %g %g %1',[11 1]);
		elseif strcmp(A,' Estimate'),
			continue;
		elseif strcmp(A,'Estimation'),
			A=fscanf(fid,'%s',1); %Read Variance
			A=fscanf(fid,'%g %g',[params(1) params(2)*params(3)]);
			B=A(1,:); B=reshape(B,[params(3),params(2)])';
			E=A(2,:); E=reshape(E,[params(3),params(2)])';
		else
			%do nothing
		end
	end
	fclose(fid);
end
