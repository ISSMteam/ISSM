function md=postqmu(md)
%POSTQMU - Deal with Dakota output results in files
%
%   Usage:
%      md = postqmu(md)

%  check to see if dakota returned errors in the err file
qmuerrfile=[md.miscellaneous.name '.qmu.err'];

if exist(qmuerrfile,'file')
   fide=fopen(qmuerrfile,'r');
   fline=fgetl(fide);
   if ischar(fline)
       while ischar(fline)
           disp(sprintf('%s',fline));
           fline=fgetl(fide);
       end
       warning(['Dakota returned error in ''' qmuerrfile ''' file. QMU directory retained.'])
   end
   status=fclose(fide);
end

%parse inputs and results from dakota
qmuinfile=[md.miscellaneous.name '.qmu.in'];
qmuoutfile=[md.miscellaneous.name '.qmu.out'];
[method,dresp_out,scm,pcm,srcm,prcm]=dakota_out_parse(qmuoutfile);
dakotaresults.dresp_out=dresp_out;
dakotaresults.scm      =scm;
dakotaresults.pcm      =pcm;
dakotaresults.srcm     =srcm;
dakotaresults.prcm     =prcm;

if exist('dakota_tabular.dat','file')
    [method,dresp_dat                  ]=dakota_out_parse('dakota_tabular.dat');
    dakotaresults.dresp_dat=dresp_dat;
end

if md.qmu.output & strcmpi(md.qmu.statistics.method(1).name,'None'),
	if strcmpi(md.qmu.method.method,'nond_sampling'),
		dakotaresults.modelresults={};
		md2=md;
		md2.qmu.isdakota=0;
		for i=1:md2.qmu.method.params.samples,
			disp(['Reading qmu file ' md2.miscellaneous.name '.outbin.' num2str(i)]);
			md2=loadresultsfromdisk(md2,[md2.miscellaneous.name '.outbin.' num2str(i)]);
			dakotaresults.modelresults{end+1}=md2.results;
		end
	end
end
if ~strcmpi(md.qmu.statistics.method(1).name,'None'),
	md.qmu.isdakota=0;
	md=loadresultsfromdisk(md,[md.miscellaneous.name,'.stats']);
	md.qmu.isdakota=1;
end

%put dakotaresults in their right location.
md.results.dakota=dakotaresults;

%  move all the individual function evalutations into zip files
if ~md.qmu.isdakota,
	system('zip -mq params.in.zip params.in.[1-9]*');
	system('zip -mq results.out.zip results.out.[1-9]*');
	system('zip -mq matlab.out.zip matlab*.out.[1-9]*');
end
