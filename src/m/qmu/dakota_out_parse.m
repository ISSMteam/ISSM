function [method,dresp,scm,pcm,srcm,prcm]=dakota_out_parse(filei) % {{{
%DAKOTA_OUT_PARSE - Read a Dakota .out or .dat output file and parse it.
%
%   Usage:
%      [method,dresp,scm,pcm,srcm,prcm]=dakota_out_parse(filei)
%
%   where the required input is,
%      filei         (character, name of .out file)
%
%   the required output is,
%      method        (character, dakota method name)
%      dresp         (structure array, responses)
%
%   and the optional output is,
%      scm           (double array, simple correlation matrix)
%      pcm           (double array, partial correlation matrix)
%      srcm          (double array, simple rank correlation matrix)
%      prcm          (double array, partial rank correlation matrix)
%
%   The filei will be prompted for if empty. The fields of dresp are particular 
%   to the data contained within the file. The scm, pcm, srcm, and prcm are 
%   output by dakota only for the sampling methods.
%
%   This function reads a dakota .out output file and parses it into the MATLAB 
%   workspace. It operates in a content-driven fashion: it skips the 
%   intermediate data and then parses whatever output data it encounters in the 
%   order in which it exists in the file, rather than searching for data based 
%   on the particular method (this makes it independent of method). It also can 
%   read and parse the .dat tabular_output file.
%
%   This data would typically be used for plotting and other postprocessing 
%   within MATLAB or Excel.
%
%   "Copyright 2009, by the California Institute of Technology. ALL RIGHTS 
%   RESERVED. United States Government Sponsorship acknowledged. Any commercial 
%   use must be negotiated with the Office of Technology Transfer at the 
%   California Institute of Technology. (NTR 47078)
%
%   This software may be subject to U.S. export control laws. By accepting this 
%   software, the user agrees to comply with all applicable U.S. export laws 
%   and regulations. User has the responsibility to obtain export licenses, or 
%   other export authority as may be required before exporting such information 
%   to foreign countries or providing access to foreign persons."
%
	if ~nargin
		help dakota_out_parse
		return
	end

	if ~exist('filei' ,'var') || isempty(filei)
		filei=input('Input file?  ','s');
	end
	fidi=fopen(sprintf('%s',filei),'r');
	if (fidi < 0)
		error('''%s'' could not be opened',filei);
	end

	%%  check the first line for the Dakota tabular output file

	method=[];
	fline=fgetl(fidi);
	if ~ischar(fline)
		error('File ''%s'' is empty',filei);
	end

	if strncmpi(fline,'%eval_id',8)
		method='unknown';
		[dresp]=dak_tab_out(fidi,fline);
		return
	else
		fseek(fidi,0,'bof');
	end

	%%  loop through the file to find the Dakota method name

	[fline]=findline(fidi,'method');
	if ~ischar(fline)
		%do nothing
	else
		% display(['  ' deblank(fline)]);
		[ntokens,tokens]=fltokens(fline);
		%dakota version >6
		if strcmp(tokens{1}{1}(7),',')
			fline=fgetl(fidi);
			[ntokens,tokens]=fltokens(fline);
			method=tokens{1}{1};
			display(sprintf('Dakota method = ''%s''',method));
		elseif strcmp(tokens{1}{1}(7),'N')
			[fline]=findline(fidi,'methodName = ');
			[ntokens,tokens]=fltokens(fline);
			method=tokens{1}{3};
			display(sprintf('Dakota methodName = ''%s''',method));
		end
	end

	dresp=struct([]);
	scm =struct([]);
	pcm =struct([]);
	srcm=struct([]);
	prcm=struct([]);

	%%  loop through the file to find the function evaluation summary

	fline='';
	[nfeval]=nfeval_read(fidi,fline);
	fline=fgetl(fidi);

	%%  process each results section based on content of the file

	while ischar(fline)
		%     ipos=ftell(fidi);
		if     isempty(fline)
		elseif strncmp(fline,'<<<<< Function evaluation summary',33)
			[nfeval]=nfeval_read(fidi,fline);
		elseif strncmp(fline,'Statistics based on ',20)
			[nsamp]=nsamp_read(fidi,fline);
		elseif strncmp(fline,'Moments for each response function',34)
			[dresp]=moments_read(fidi,dresp,fline);
		elseif strncmp(fline,'Moment-based statistics for each response function',50)
			[dresp]=mbstats_read(fidi,dresp,fline);
		elseif strncmp(fline,'95% confidence intervals for each response function',51)
			[dresp]=cis_read(fidi,dresp,fline);
		elseif strncmp(fline,'Probabilities for each response function',40) || ...
				strncmp(fline,'Level mappings for each response function',41)
			[dresp]=cdfs_read(fidi,dresp,fline);
		elseif strncmp(fline,'Probability Density Function (PDF) histograms for each response function',72)
			[dresp]=pdfs_read(fidi,dresp,fline);
		elseif strncmp(fline,'Simple Correlation Matrix',25)
			[scm]=corrmat_read(fidi,'Simple Correlation Matrix',fline);
		elseif strncmp(fline,'Partial Correlation Matrix',26)
			[pcm]=corrmat_read(fidi,'Partial Correlation Matrix',fline);
		elseif strncmp(fline,'Simple Rank Correlation Matrix',30)
			[srcm]=corrmat_read(fidi,'Simple Rank Correlation Matrix',fline);
		elseif strncmp(fline,'Partial Rank Correlation Matrix',31)
			[prcm]=corrmat_read(fidi,'Partial Rank Correlation Matrix',fline);
		elseif strncmp(fline,'MV Statistics for ',18)
			[dresp]=mvstats_read(fidi,dresp,fline);
		elseif strncmp(fline,'<<<<< Best ',11)
			[dresp]=best_read(fidi,dresp,fline);
		elseif strncmp(fline,'The following lists volumetric uniformity measures',50)
			[dresp]=vum_read(fidi,dresp,fline);
		elseif strncmp(fline,'<<<<< Iterator ',15) && ...
				(length(fline) > 26) && ...
				~isempty(strfind(fline(16:end),' completed.'))
			[method]=itcomp_read(fidi,fline);
		elseif strncmp(fline,'-----',5)
		elseif strncmp(fline,'<<<<< Environment execution completed',37)
			break;
		else
			display(['Unexpected line: ' deblank(fline)]);
		end
		fline=fgetl(fidi);
		%     fseek(fidi,ipos,'bof');
	end

	%%  loop through the file to verify the end

	% [fline]=findline(fidi,'<<<<< Single Method Strategy completed');
	% if ~ischar(fline)
	%     return
	% end
	display('End of file successfully reached');
	fclose(fidi);

end % }}}
function [dresp]=dak_tab_out(fidi,fline) % {{{
%%  function to parse the dakota tabular output file

	display('Reading Dakota tabular output file');

	%  process column headings of matrix (skipping eval_id)

	[ntokens,tokens]=fltokens(fline);

	if strncmpi(fline,'%eval_id interface',18) % Dakota versions >= 6
		offset=2;
	else % Dakota versions < 6
		offset=1;
	end
	desc=cell(1,ntokens-offset);
	data=zeros(1,ntokens-offset);

	for i=1:ntokens-offset
		desc(1,i)=cellstr(tokens{1}{i+offset});
	end
	display(sprintf('Number of columns (Dakota V + R) = %d',ntokens-2));

	%  process rows of matrix

	nrow=0;
	while 1
		fline=fgetl(fidi);

		if ~ischar(fline) || isempty(fline)
			break;
		end
		[ntokens,tokens]=fltokens(fline);

		%  add row values to matrix (skipping eval_id)

		nrow=nrow+1;
		for i=1:ntokens-offset
			data(nrow,i)=tokens{1}{i+offset};
		end
	end
	display(sprintf('Number of rows (Dakota func evals) = %d',nrow));

	%  calculate statistics

	%  since normfit doesn't have a dim argument, and matlab isvector is true
	%  for a 1xn matrix, handle the case of one row explicitly
	if (size(data,1) > 1)
		%dmean  =mean   (data);
		%dstddev=std    (data,0);
		[dmean,dstddev,dmeanci,dstddevci]=...
			normfit_issm(data,0.05);
	else
		dmean    =zeros(1,size(data,2));
		dstddev  =zeros(1,size(data,2));
		dmeanci  =zeros(2,size(data,2));
		dstddevci=zeros(2,size(data,2));
		for i=1:size(data,2)
			[dmean(1,i),dstddev(1,i),dmeanci(:,i),dstddevci(:,i)]=...
				normfit_issm(data(:,i),0.05);
		end
	end

	dmin   =min(data,[],1);
	dquart1=prctile_issm(data,25,1);
	dmedian=median(data,1);
	dquart3=prctile_issm(data,75,1);
	dmax   =max(data,[],1);
	dmin95 =prctile_issm(data,5,1);
	dmax95 =prctile_issm(data,95,1);

	%  same as Dakota scm, Excel correl
	dcorrel=corrcoef(data);

	%  divide the data into structures for consistency

	for i=1:length(desc)
		dresp(i).descriptor=char(desc(i));
		dresp(i).sample    =data(:,i);
		dresp(i).mean      =dmean(i);
		dresp(i).stddev    =dstddev(i);
		dresp(i).meanci    =dmeanci(:,i);
		dresp(i).stddevci  =dstddevci(:,i);
		dresp(i).min       =dmin(i);
		dresp(i).quart1    =dquart1(i);
		dresp(i).median    =dmedian(i);
		dresp(i).quart3    =dquart3(i);
		dresp(i).max       =dmax(i);
		dresp(i).dmin95    =dmin95(i);
		dresp(i).dmax95    =dmax95(i);
	end

	%  draw box plot

	% figure
	% subplot(2,1,1)
	% plot_boxplot(dresp);

	%  draw normal probability plot

	% subplot(2,1,2)
	% plot_normplot(dresp);

end % }}}
function [nfeval]=nfeval_read(fidi,fline) % {{{
%%  function to find and read the number of function evaluations

	if ~exist('fline','var') || isempty(fline) || ~ischar(fline)
		[fline]=findline(fidi,'<<<<< Function evaluation summary');
		if ~ischar(fline)
			nfeval = 0;
			return
		end
	end

	[ntokens,tokens]=fltokens(fline);
	nfeval=tokens{1}{5};
	display(sprintf('  Dakota function evaluations = %d',nfeval));

end % }}}
function [nsamp]=nsamp_read(fidi,fline) % {{{
%%  function to find and read the number of samples

	if ~exist('fline','var') || isempty(fline) || ~ischar(fline)
		[fline]=findline(fidi,'Statistics based on ');
		if ~ischar(fline)
			return
		end
	end

	[ntokens,tokens]=fltokens(fline);
	nsamp=tokens{1}{4};
	display(sprintf('  Dakota samples = %d',nsamp));

end % }}}
function [dresp]=moments_read(fidi,dresp,fline) % {{{
%%  function to find and read the moments

	if ~exist('fline','var') || isempty(fline) || ~ischar(fline)
		[fline]=findline(fidi,'Moments for each response function');
		if ~ischar(fline)
			return
		end
	end

	display('Reading moments for response functions:');

	while 1
		fline=fgetl(fidi);
		if isempty(fline)
			break;
		end
		[ntokens,tokens]=fltokens(fline);

		%  add new response function and moments

		dresp(end+1).descriptor=tokens{1}{ 1};
		display(sprintf('  %s',dresp(end).descriptor));
		dresp(end  ).mean      =tokens{1}{ 4};
		dresp(end  ).stddev    =tokens{1}{ 8};
		dresp(end  ).coefvar   =tokens{1}{13};
	end

	display(sprintf('  Number of Dakota response functions = %d',...
		length(dresp)));

end % }}}
function [dresp]=mbstats_read(fidi,dresp,fline) % {{{
%%  function to find and read the moment-based statistics

	if ~exist('fline','var') || isempty(fline) || ~ischar(fline)
		[fline]=findline(fidi,'Moment-based statistics for each response function');
		if ~ischar(fline)
			return
		end
	end

	display('Reading moment-based statistics for response functions:');

	%  skip column headings of moment-based statistics

	fline=fgetl(fidi);

	while 1
		fline=fgetl(fidi);
		if isempty(fline)
			break;
		end
		[ntokens,tokens]=fltokens(fline);

		%  add new response function and moment-based statistics

		dresp(end+1).descriptor=tokens{1}{ 1};
		display(sprintf('  %s',dresp(end).descriptor));
		dresp(end  ).mean      =tokens{1}{ 2};
		dresp(end  ).stddev    =tokens{1}{ 3};
		dresp(end  ).skewness  =tokens{1}{ 4};
		dresp(end  ).kurtosis  =tokens{1}{ 5};
	end

	display(sprintf('  Number of Dakota response functions = %d',...
		length(dresp)));

end % }}}
function [dresp]=cis_read(fidi,dresp,fline) % {{{
%%  function to find and read the confidence intervals

	if ~exist('fline','var') || isempty(fline) || ~ischar(fline)
		[fline]=findline(fidi,...
			'95% confidence intervals for each response function');
		if ~ischar(fline)
			return
		end
	end

	display('Reading 95% confidence intervals for response functions:');

	while 1
		fline=fgetl(fidi);
		if isempty(fline)
			break;
		end
		[ntokens,tokens]=fltokens(fline);

		%  check for column headings in Dakota 5.2

		if (ntokens == 4)
			fline=fgetl(fidi);
			if isempty(fline)
				break;
			end
			[ntokens,tokens]=fltokens(fline);
		end

		%  find response function associated with confidence intervals

		idresp=0;
		for i=1:length(dresp)
			if strcmpi(tokens{1}{ 1},dresp(i).descriptor)
				idresp=i;
				break;
			end
		end
		if ~idresp
			idresp=length(dresp)+1;
			dresp(idresp).descriptor=tokens{1}{ 1};
			display(sprintf('  %s',dresp(idresp).descriptor));
		end

		%  add confidence intervals to response functions

		if (ntokens == 14)
			dresp(i).meanci  (1,1)=tokens{1}{ 5};
			dresp(i).meanci  (2,1)=tokens{1}{ 6};
			dresp(i).stddevci(1,1)=tokens{1}{12};
			dresp(i).stddevci(2,1)=tokens{1}{13};
		else
			dresp(i).meanci  (1,1)=tokens{1}{ 2};
			dresp(i).meanci  (2,1)=tokens{1}{ 3};
			dresp(i).stddevci(1,1)=tokens{1}{ 4};
			dresp(i).stddevci(2,1)=tokens{1}{ 5};
		end
	end

	display(sprintf('  Number of Dakota response functions = %d',...
		length(dresp)));

end % }}}
function [dresp]=cdfs_read(fidi,dresp,fline) % {{{
%%  function to find and read the cdf's

	if ~exist('fline','var') || isempty(fline) || ~ischar(fline)
		[fline]=findline(fidi,'Probabilities for each response function');
		if ~ischar(fline)
			[fline]=findline(fidi,'Level mappings for each response function');
			if ~ischar(fline)
				return
			end
		end
	end

	display('Reading CDF''s for response functions:');

	while ischar(fline) && ~isempty(fline)
		fline=fgetl(fidi);
		if ~ischar(fline)
			break;
		end

		%  process header line of cdf

		while ischar(fline) && ~isempty(fline)
			[ntokens,tokens]=fltokens(fline);

			%  find response function associated with cdf

			idresp=0;
			for i=1:length(dresp)
				if strcmpi(tokens{1}{ 6},dresp(i).descriptor)
					idresp=i;
					break;
				end
			end
			if ~idresp
				idresp=length(dresp)+1;
				dresp(idresp).descriptor=tokens{1}{ 6};
				display(sprintf('  %s',dresp(idresp).descriptor));
			end

			%  skip column headings of cdf

			fline=fgetl(fidi);
			fline=fgetl(fidi);

			%  read and add cdf table to response function

			fline=fgetl(fidi);
			icdf=0;
			while ischar(fline) && ~isempty(fline) && ...
					~strncmpi(fline,'Cumulative Distribution Function',32)
				[ntokens,tokens]=fltokens(fline);
				icdf=icdf+1;
				dresp(idresp).cdf(icdf,1:4)=NaN;
				%  in later versions of Dakota, uncalculated columns are now blank
				itoken=0;
				for i=1:length(fline)/19
					if ~isempty(deblank(fline((i-1)*19+1:i*19)))
						itoken=itoken+1;
						dresp(idresp).cdf(icdf,i)=tokens{1}{itoken};
					end
				end
				fline=fgetl(fidi);
			end
		end
	end

	display(sprintf('  Number of Dakota response functions = %d',...
		length(dresp)));

end % }}}
function [dresp]=pdfs_read(fidi,dresp,fline) % {{{
%%  function to find and read the pdf's

	if ~exist('fline','var') || isempty(fline) || ~ischar(fline)
		[fline]=findline(fidi,'Probability Density Function (PDF) histograms for each response function');
		if ~ischar(fline)
			return
		end
	end

	display('Reading PDF''s for response functions:');

	while ischar(fline) && ~isempty(fline)
		fline=fgetl(fidi);
		if ~ischar(fline)
			break;
		end

		%  process header line of pdf

		while ischar(fline) && ~isempty(fline)
			[ntokens,tokens]=fltokens(fline);

			%  find response function associated with pdf

			idresp=0;
			for i=1:length(dresp)
				if strcmpi(tokens{1}{ 3},dresp(i).descriptor)
					idresp=i;
					break;
				end
			end
			if ~idresp
				idresp=length(dresp)+1;
				dresp(idresp).descriptor=tokens{1}{ 3};
				display(sprintf('  %s',dresp(idresp).descriptor));
			end

			%  skip column headings of pdf

			fline=fgetl(fidi);
			fline=fgetl(fidi);

			%  read and add pdf table to response function

			fline=fgetl(fidi);
			ipdf=0;
			while ischar(fline) && ~isempty(fline) && ...
					~strncmpi(fline,'PDF for', 7)
				[ntokens,tokens]=fltokens(fline);
				ipdf=ipdf+1;
				dresp(idresp).pdf(ipdf,1:3)=NaN;
				for i=1:3
					dresp(idresp).pdf(ipdf,i)=tokens{1}{i};
				end
				fline=fgetl(fidi);
			end
		end
	end

	display(sprintf('  Number of Dakota response functions = %d',...
		length(dresp)));

end % }}}
function [cmat]=corrmat_read(fidi,cmstr,fline) % {{{
%%  function to find and read a correlation matrix

	if ~exist('fline','var') || isempty(fline) || ~ischar(fline)
		[fline]=findline(fidi,cmstr);
		if ~ischar(fline)
			cmat=struct([]);
			return
		end
	end

	display(['Reading ''' fline '''']);

	cmat.title=fline;

	while ~isempty(fline)
		fline=fgetl(fidi);
		if ~ischar(fline)
			break;
		end

		%  process column headings of matrix

		[ntokens,tokens]=fltokens(fline);
		cmat.column=cell(1,ntokens);
		cmat.row   =cell(1,1);
		cmat.matrix=zeros(1,ntokens);

		for i=1:ntokens
			cmat.column(1,i)=cellstr(tokens{1}{i});
		end

		%  process rows of matrix, reading until blank line

		nrow=0;
		while 1
			fline=fgetl(fidi);
			if isempty(fline)
				break;
			end
			[ntokens,tokens]=fltokens(fline);

			%  add row heading to matrix

			nrow=nrow+1;
			cmat.row   (nrow,1)=cellstr(tokens{1}{1});

			%  add row values to matrix

			for i=2:ntokens
				cmat.matrix(nrow,i-1)=tokens{1}{i};
			end
		end
	end

end % }}}
function [dresp]=mvstats_read(fidi,dresp,fline) % {{{
%%  function to find and read the MV statistics

	if ~exist('fline','var') || isempty(fline) || ~ischar(fline)
		[fline]=findline(fidi,'MV Statistics for ');
		if ~ischar(fline)
			return
		end
	end

	display('Reading MV statistics for response functions:');
	ndresp=0;

	while ischar(fline) && ~isempty(fline) && ...
			strncmpi(fline,'MV Statistics for ',18)

		%  add new response function and moments

		[ntokens,tokens]=fltokens(fline);
		dresp(end+1).descriptor=tokens{1}{4};
		display(sprintf('  %s',dresp(end).descriptor));
		fline=fgetl(fidi);
		[ntokens,tokens]=fltokens(fline);
		dresp(end  ).mean      =tokens{1}{5};
		fline=fgetl(fidi);
		[ntokens,tokens]=fltokens(fline);
		dresp(end  ).stddev    =tokens{1}{7};

		%  read and add importance factors to response function

		idvar=0;
		fline=fgetl(fidi);
		if ~ischar(fline)
			break;
		end

		while ischar(fline) && ~isempty(fline) && ...
				strncmpi(fline,'  Importance Factor for variable ',33)
			[ntokens,tokens]=fltokens(fline);
			idvar=idvar+1;
			dresp(end).var   (idvar,1)=cellstr(tokens{1}{ 5});
			dresp(end).impfac(idvar,1)=        tokens{1}{ 7};
			if (ntokens >= 10)
				dresp(end).sens  (idvar,1)=        tokens{1}{10};
			else
				dresp(end).sens  (idvar,1)=NaN;
			end

			fline=fgetl(fidi);
		end

		%  if importance factors missing, skip to cdf

		if ~idvar
			display('    Importance Factors not available');
			dresp(end).var   ={};
			dresp(end).impfac=[];
			dresp(end).sens  =[];
			while ischar(fline) && ...
					~strncmpi(fline,'Cumulative Distribution Function',32) && ...
					~strncmpi(fline,'MV Statistics for ',18) && ...
					~strncmp (fline,'-',1)
				fline=fgetl(fidi);
			end
		end

		%  process header line of cdf

		icdf=0;

		while ischar(fline) && ~isempty(fline) && ...
				strncmpi(fline,'Cumulative Distribution Function',32)
			[ntokens,tokens]=fltokens(fline);

			%  find response function associated with cdf

			idresp=0;
			for i=1:length(dresp)
				if strcmpi(tokens{1}{ 6},dresp(i).descriptor)
					idresp=i;
					break;
				end
			end
			if ~idresp
				idresp=length(dresp)+1;
				dresp(idresp).descriptor=tokens{1}{ 6};
				display(sprintf('  %s',dresp(idresp).descriptor));
			end

			%  skip column headings of cdf

			fline=fgetl(fidi);
			fline=fgetl(fidi);

			%  read and add cdf table to response function

			fline=fgetl(fidi);
			while ~isempty(fline) && ...
					~strncmpi(fline,'MV Statistics for ',18) && ...
					~strncmp (fline,'-',1)
				[ntokens,tokens]=fltokens(fline);
				icdf=icdf+1;
				dresp(idresp).cdf(icdf,1)=tokens{1}{1};
				dresp(idresp).cdf(icdf,2)=tokens{1}{2};
				if (ntokens == 4)
					dresp(idresp).cdf(icdf,3)=tokens{1}{3};
					dresp(idresp).cdf(icdf,4)=tokens{1}{4};
				else
					dresp(idresp).cdf(icdf,3)=NaN;
					dresp(idresp).cdf(icdf,4)=NaN;
				end
				fline=fgetl(fidi);
			end
		end

		%  if cdf missing, skip to end of response function

		if ~icdf
			display('    Cumulative Distribution Function not available');
			dresp(ndresp).cdf=[];
			while ischar(fline) && ...
					~strncmpi(fline,'MV Statistics for ',18) && ...
					~strncmp (fline,'-',1)
				fline=fgetl(fidi);
			end
		end

	end

	display(sprintf('  Number of Dakota response functions = %d',...
		length(dresp)));

end % }}}
function [dresp]=best_read(fidi,dresp,fline) % {{{
%%  function to find and read the best evaluation

	if ~exist('fline','var') || isempty(fline) || ~ischar(fline)
		[fline]=findline(fidi,'<<<<< Best ');
		if ~ischar(fline)
			return
		end
	end

	if isempty(dresp)
		dresp(end+1).best=[];
	end
	display('Reading values for best function evaluation:');

	while ischar(fline) && ~isempty(fline) && ...
			strncmpi(fline,'<<<<< Best ',11)
		[ntokens,tokens]=fltokens(fline);

		%  read and add best parameter(s)

		if     strncmpi(cellstr(tokens{1}{3}),'parameter', 9)
			display(['  ' deblank(fline)]);

			fline=fgetl(fidi);
			dresp.best.param     =[];
			dresp.best.descriptor={};

			while ischar(fline) && ~isempty(fline) && ...
					~strncmpi(fline,'<<<<< Best ',11)
				[ntokens,tokens]=fltokens(fline);
				dresp.best.param     (end+1,1)=        tokens{1}{1};
				%dresp.best.descriptor(end+1,1)=cellstr(tokens{1}{2}); this line is wrong, there is no descriptor
				fline=fgetl(fidi);
			end

			%  read and add best objective function(s)

		elseif strncmpi(cellstr(tokens{1}{3}),'objective', 9) && ...
				strncmpi(cellstr(tokens{1}{4}),'function' , 8)
			display(['  ' deblank(fline)]);

			fline=fgetl(fidi);
			dresp.best.of=[];

			while ischar(fline) && ~isempty(fline) && ...
					~strncmpi(fline,'<<<<< Best ',11)
				[ntokens,tokens]=fltokens(fline);
				dresp.best.of(end+1,1)=        tokens{1}{1};
				fline=fgetl(fidi);
			end

			%  read and add best residual norms

		elseif strncmpi(cellstr(tokens{1}{3}),'residual', 8) && ...
				strncmpi(cellstr(tokens{1}{4}),'norm'    , 4)
			display(['  ' deblank(fline)]);
			dresp.best.norm   =        tokens{1}{ 6};
			dresp.best.hnormsq=        tokens{1}{11};

			fline=fgetl(fidi);

			while ischar(fline) && ~isempty(fline) && ...
					~strncmpi(fline,'<<<<< Best ',11)
				fline=fgetl(fidi);
			end

			%  read and add best residual term(s)

		elseif strncmpi(cellstr(tokens{1}{3}),'residual', 8) && ...
				strncmpi(cellstr(tokens{1}{4}),'term'    , 4)
			display(['  ' deblank(fline)]);

			fline=fgetl(fidi);
			dresp.best.res=[];

			while ischar(fline) && ~isempty(fline) && ...
					~strncmpi(fline,'<<<<< Best ',11)
				[ntokens,tokens]=fltokens(fline);
				dresp.best.res(end+1,1)=        tokens{1}{1};
				fline=fgetl(fidi);
			end

			%  read and add best constraint value(s)

		elseif strncmpi(cellstr(tokens{1}{3}),'constraint',10) && ...
				strncmpi(cellstr(tokens{1}{4}),'value'     , 5)
			display(['  ' deblank(fline)]);

			fline=fgetl(fidi);
			dresp.best.nc=[];

			while ischar(fline) && ~isempty(fline) && ...
					~strncmpi(fline,'<<<<< Best ',11)
				[ntokens,tokens]=fltokens(fline);
				dresp.best.nc(end+1,1)=        tokens{1}{1};
				fline=fgetl(fidi);
			end

			%  read and add best data captured

		elseif strncmpi(cellstr(tokens{1}{3}),'data'    , 4) && ...
				strncmpi(cellstr(tokens{1}{4}),'captured', 8)
			display(['  ' deblank(fline)]);
			dresp.best.eval=        tokens{1}{8};

			fline=fgetl(fidi);

			while ischar(fline) && ~isempty(fline) && ...
					~strncmpi(fline,'<<<<< Best ',11)
				fline=fgetl(fidi);
			end

			%  read until next best or blank or end

		else
			display(['  ' deblank(fline) '  (ignored)']);

			fline=fgetl(fidi);

			while ischar(fline) && ~isempty(fline) && ...
					~strncmpi(fline,'<<<<< Best ',11)
				fline=fgetl(fidi);
			end
		end
	end

end % }}}
function [dresp]=vum_read(fidi,dresp,fline) % {{{
%%  function to find and read the volumetric uniformity measures

	if ~exist('fline','var') || isempty(fline) || ~ischar(fline)
		[fline]=findline(fidi,'The following lists volumetric uniformity measures');
		if ~ischar(fline)
			return
		end
	end

	if isempty(dresp)
		dresp(end+1).vum=[];
	end
	display('Reading measures for volumetric uniformity');

	fline=fgetl(fidi);
	fline=fgetl(fidi);

	while ischar(fline) && ~isempty(fline)
		[ntokens,tokens]=fltokens(fline);
		switch lower(tokens{1}{1})
			case 'chi'
				dresp.vum.chi=tokens{1}{4};
			case 'd'
				dresp.vum.d  =tokens{1}{4};
			case 'h'
				dresp.vum.h  =tokens{1}{4};
			case 'tau'
				dresp.vum.tau=tokens{1}{4};
		end
		fline=fgetl(fidi);
	end

end % }}}
function [method]=itcomp_read(fidi,fline) % {{{
%%  function to find and read the iterator completion

	if ~exist('fline','var') || isempty(fline) || ~ischar(fline)
		while 1
			[fline]=findline(fidi,'<<<<< Iterator ');
			if ~ischar(fline)
				return
			end
			if (length(fline) > 26) && ...
					~isempty(strfind(fline(16:end),' completed.'))
				break
			end
		end
	end

	[ntokens,tokens]=fltokens(fline);
	method=tokens{1}{3};
	display(sprintf('Dakota iterator ''%s'' completed',method));

end % }}}
function [fline]=findline(fidi,string) % {{{
%%  function to find a file line starting with a specified string

	ipos=ftell(fidi);

	while 1
		fline=fgetl(fidi);
		if ~ischar(fline)
			break;
		else
			if (strncmpi(fline,string,length(string)))
				return;
			end
		end
	end

	%  issue warning and reset file position

	warning('findline:str_not_found',...
		'String ''%s'' not found in file',string);
	fseek(fidi,ipos,'bof');

end % }}}
function [ntokens,tokens]=fltokens(fline) % {{{
%%  function to parse a file line into tokens

	if ~ischar(fline)
		ntokens=-1;
		tokens={};
		return;
	end
	if isempty(fline)
		ntokens=0;
		tokens={};
		return;
	end

	strings=textscan(fline,'%s','delimiter',' :');
	%for i=1:length(strings{1})
	%    display(sprintf('i=%d; strings{1}{%d}=%s',i,i,strings{1}{i}))
	%end
	ntokens=0;
	tokens{1}{length(strings)}='';

	for i=1:length(strings{1})
		if isempty(strings{1}{i})
			continue
		end
		ntokens=ntokens+1;
		inum=sscanf(strings{1}{i},'%f');
		if isempty(inum)
			tokens{1}{ntokens}=strings{1}{i};
			%        display(sprintf('i=%d; tokens{1}{%d}=%s',...
			%            i,ntokens,tokens{1}{ntokens}))
		else
			tokens{1}{ntokens}=inum;
			%        display(sprintf('i=%d; tokens{1}{%d}=%f',...
			%            i,ntokens,tokens{1}{ntokens}))
		end
	end

end % }}}
