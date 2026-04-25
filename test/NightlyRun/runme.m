function varargout=runme(varargin)
%RUNME - test deck for ISSM nightly runs
%
%   In a test deck directory (for example, test/NightlyRun) the following
%   command will launch all existing tests,
%
%      >> runme
%
%   To run tests 101 and 102,
%
%      >> runme('id',[101 102])
%
%   Available options:
%      'id'            followed by the list of test ID's or names to run
%      'exclude'       followed by the list of test ID's or names to exclude
%      'benchmark'     'all'         : (all of them)
%                      'nightly'     : (nightly run
%                      'validation'  : (validation)
%                      'adolc'       : validation of adolc tests
%                      'eismint'     : validation of eismint tests
%                      'ismip'       : validation of ismip-hom tests
%                      'mesh'        : validation of mesh tests
%                      'qmu'         : validation of dakota tests
%                      'referential' : validation of referential tests
%                      'slc'         : validation of slc tests
%                      'thermal'     : validation of thermal tests
%                      'tranforcing' : validation of transient forcing tests
%      'procedure'     'check' :   run the test (default)
%                      'update':   update the archive
%                      'valgrind': check for memory leaks (default value of md.debug.valgrind needs to be changed manually)
%                      'ncExport': export netCDF file
%      'stoponerror'   1 or 0
%      'quitonerror'   1 or 0
%
%   Usage:
%      runme(varargin);
%
%   Examples:
%      runme;
%      runme('exclude',101);
%      runme('id', 102, 'procedure', 'update');
%      runme('procedure', 'valgrind', 'stoponerror', 1, 'exclude', IdFromString('Dak'));
%
%   NOTE:
%   - Will only run test scripts whose names explicitly follow the convention,
%
%         test<id>.m
%
%   where <id> is any integer.
%

%Check inputs
% {{{
if nargout>1
	help runme
	error('runme error message: bad usage');
end

%recover options
options=pairoptions(varargin{:});
% }}}

%Process options
%Get benchmark {{{
benchmark=getfieldvalue(options,'benchmark','nightly');
if ~ismember(benchmark,{'all','nightly','ismip','eismint','thermal','mesh','validation','tranforcing','adolc','slc','qmu'})
	disp('runme warning: benchmark not supported, defaulting to test ''nightly''')
	benchmark='nightly';
end
% }}}
%Get procedure {{{
procedure=getfieldvalue(options,'procedure','check');
if ~ismember(procedure,{'check','update','valgrind','ncExport'})
	disp('runme warning: procedure not supported, defaulting to test ''check''')
	procedure='check';
end
% }}}
%Get output {{{
output=getfieldvalue(options,'output','none');
if ~ismember(output,{'nightly','none'})
	disp('runme warning: output not supported, defaulting to test ''none''')
	output='none';
end
% }}}
%Get rank and numprocs for multi-threaded runs  {{{
rank=getfieldvalue(options,'rank',1);
numprocs=getfieldvalue(options,'numprocs',1);
if (numprocs<rank), numprocs=1; end
% }}}
%Get IDs (create a list of all the test files in this directory that match a certain naming scheme) {{{
flist=dir; %use dir, as it seems to act OS independent
list_ids=[];
for i=1:numel(flist),
	fname=flist(i).name;
	if (any(fname=='.')), %before split, check that file name contains '.'
		ftokens=strsplit(fname,'.'); %tokenize file name on '.'
		if (regexp(ftokens{1},'^test[0-9]+$') &... %basename must start with 'test' and end with an integer
			strcmp(ftokens{end},'m') ... %extension (less '.') must be 'm'
		),
			id=sscanf(ftokens{1},'test%d');
			if isempty(id),
				disp(['WARNING: ignore file ' flist(i).name]);
			else
				list_ids(end+1)=id; %keep test id only (strip 'test' and '.m')
			end
		end
	end
end
[i1,i2]=parallelrange(rank,numprocs,length(list_ids)); %get tests for this CPU only
list_ids=list_ids(i1:i2);

test_ids=getfieldvalue(options,'id',list_ids);
test_ids=intersect(test_ids,list_ids);
% }}}
%Get excluded tests {{{
exclude_ids=getfieldvalue(options,'exclude',[]);
exclude_ids=[exclude_ids];
pos=find(ismember(test_ids,exclude_ids));
test_ids(pos)=[];
% }}}
%Process IDs according to benchmarks{{{
if strcmpi(benchmark,'nightly'),
	test_ids=intersect(test_ids,[1:999]);
elseif strcmpi(benchmark,'validation'),
	test_ids=intersect(test_ids,[1001:1999]);
elseif strcmpi(benchmark,'ismip'),
	test_ids=intersect(test_ids,[1101:1199]);
elseif strcmpi(benchmark,'eismint'),
	test_ids=intersect(test_ids,[1201:1299]);
elseif strcmpi(benchmark,'thermal'),
	test_ids=intersect(test_ids,[1301:1399]);
elseif strcmpi(benchmark,'mesh'),
	test_ids=intersect(test_ids,[1401:1499]);
elseif strcmpi(benchmark,'tranforcing'),
	test_ids=intersect(test_ids,[1501:1502]);
elseif strcmpi(benchmark,'referential'),
	test_ids=intersect(test_ids,[1601:1602]);
elseif strcmpi(benchmark,'slc'),
	test_ids=intersect(test_ids,[2001:2500]);
elseif strcmpi(benchmark,'adolc'),
	test_ids=intersect(test_ids,[3001:3200]);
elseif strcmpi(benchmark,'qmu'),
	test_ids=intersect(test_ids,[218 234 235 412:414 417 418 420]);
end
% }}}

%Loop over tests and launch sequence
root=pwd;
for id=test_ids,
	disp(sprintf('%s%i%s','----------------starting:',id,'-----------------------'));
	try,
		%Execute test
		cd(root);
		id_string='N/A';
		id_string=IdToName(id);
		run(['test' num2str(id)]);

		%Update archive?
		archive_name=['Archive' num2str(id) ];
		if strcmpi(procedure,'update'),
			delete(['../Archives/' archive_name '.arch'])
			for k=1:length(field_names),
				field=field_values{k};
				archwrite(['../Archives/' archive_name '.arch'],[archive_name '_field' num2str(k)], field);
			end
			disp(sprintf(['File ./../Archives/' archive_name '.arch saved\n']));

		%Check for memory leaks?
		elseif strcmpi(procedure,'valgrind'),
			fields = fieldnames(md.results);
			for i=1:numel(fields)
				if ~isfield(md.results.(fields{i}),'errlog'),
					disp(['Skipping ' fields{i}]);
					continue;
				else
					disp(['Extracting results of ' fields{i}]);
				end
				results = md.results.(fields{i});
				errlog  = cellstr(results(1).errlog);

				%Check leaks
				lines  = strfind(errlog,'definitely lost:');
				lines  = find(~cellfun(@isempty,lines));
				leaks   = 0;
				for j=1:numel(lines)
					Line    = errlog(lines(j));
					Numbers = sscanf(Line{1},'==%i==   definitely lost: %s bytes in %i blocks',[1 Inf]);
					leaks   = leaks + str2num(strrep(char(Numbers(2:end-1)),',',''));
				end
				%Check conditional jumps
				lines  = strfind(errlog,'Conditional jump or move depends on uninitialised value');
				lines  = find(~cellfun(@isempty,lines));
				jumps   = numel(lines);
				%Check invalid read/write
				lines  = strfind(errlog,'Invalid');
				lines  = find(~cellfun(@isempty,lines));
				inval  = numel(lines);
				if leaks==0,
					disp(sprintf(['SUCCESS difference: 0 < 0 test id: %i test name: %s field: valgrind mem. leaks'],id,id_string));
				else
					disp(sprintf(['ERROR   difference: %i > 0 test id: %i test name: %s field: valgrind mem. leaks'],leaks,id,id_string));
					disp('STOP');
					return;
				end
				if jumps==0,
					disp(sprintf(['SUCCESS difference: 0 < 0 test id: %i test name: %s field: valgrind cond. jumps'],id,id_string));
				else
					disp(sprintf(['ERROR   difference: %i > 0 test id: %i test name: %s field: valgrind cond. jumps'],jumps,id,id_string));
					disp('STOP');
					return;
				end
				if inval==0,
					disp(sprintf(['SUCCESS difference: 0 < 0 test id: %i test name: %s field: valgrind invalid read/write'],id,id_string));
				else
					disp(sprintf(['ERROR   difference: %i > 0 test id: %i test name: %s field: valgrind invalid read/write'],inval,id,id_string));
					disp('STOP');
					return;
				end
			end
		%Produce nc files?
		elseif strcmpi(procedure,'ncExport'),
			export_netCDF(md, ['test' num2str(id) 'ma.nc'])

		%Else: Check test
		else,
			for k=1:length(field_names),

				try,
					%Get field and tolerance
					field=field_values{k};
					fieldname=field_names{k};
					tolerance=field_tolerances{k};

					%Compare to archive
					%Our output is in the correct order (n,1) or (1,1), so we do not need to transpose again
					archive_cell=archread(['../Archives/' archive_name '.arch'],[archive_name '_field' num2str(k)]);
					archive=archive_cell{1};
					error_diff=full(max(abs(archive(:)-field(:)))/(max(abs(archive(:)))+eps)); %disp test result
					if (error_diff>tolerance | isnan(error_diff));
						disp(sprintf(['ERROR   difference: %-7.2g > %7.2g test id: %i test name: %s field: %s'],...
							error_diff,tolerance,id,id_string,fieldname));
						if(getfieldvalue(options,'stoponerror',0)), disp('STOP'); return; end
						if(getfieldvalue(options,'quitonerror',0)), disp('STOP'); quit(1); end
					else
						disp(sprintf(['SUCCESS difference: %-7.2g < %7.2g test id: %i test name: %s field: %s'],...
							error_diff,tolerance,id,id_string,fieldname));
					end

				catch me2

					%Something went wrong, print failure message
					message=getReport(me2);
					fprintf('%s',message);
					if strcmpi(output,'nightly')
						fid=fopen([issmdir() '/nightlylog/matlaberror.log'], 'at');
						fprintf(fid,'%s',message);
						fprintf(fid,'\n------------------------------------------------------------------\n');
						fclose(fid);
						disp(sprintf(['FAILURE difference: N/A test id: %i test name: %s field: %s'],id,id_string,fieldname));
					else
						disp(sprintf(['FAILURE difference: N/A test id: %i test name: %s field: %s'],id,id_string,fieldname));
						fprintf('%s',message);
						if(getfieldvalue(options,'stoponerror',0)), disp('STOP'); return; end
						if(getfieldvalue(options,'quitonerror',0)), disp('STOP'); quit(1); end
					end
					continue;
				end
			end
		end
	catch me,

		%Something went wrong, print failure message
		message=getReport(me);
		fprintf('%s',message);
		if strcmpi(output,'nightly')
			fid=fopen([issmdir() '/nightlylog/matlaberror.log'], 'at');
			fprintf(fid,'%s',message);
			fprintf(fid,'\n------------------------------------------------------------------\n');
			fclose(fid);
			disp(sprintf(['FAILURE difference: N/A test id: %i test name: %s field: %s'],id,id_string,'N/A'));
		else
			disp(sprintf(['FAILURE difference: N/A test id: %i test name: %s field: %s'],id,id_string,'N/A'));
			rethrow(me);
			if(getfieldvalue(options,'stoponerror',0)), disp('STOP'); return; end
			if(getfieldvalue(options,'quitonerror',0)), disp('STOP'); quit(1); end
		end
	end
	disp(sprintf('%s%i%s','----------------finished:',id,'-----------------------'));
end
