%PLOTOPTIONS class definition
%
%   Usage:
%      plotoptions = plotoptions(varargin)

classdef plotoptions
    properties (SetAccess=public) 
		 % {{{
		 numberofplots = 0;
		 figurenumber  = 1;
		 list          = cell(0,0);
		 %}}}
	 end
	 methods
		 function opt=plotoptions(varargin) % {{{
			 opt=buildlist(opt,varargin{:});
		 end
		 %}}}
		 function disp(opt) % {{{
			 disp(sprintf('\n%s = \n',inputname(1)));
			 disp(sprintf('   numberofplots: %i',opt.numberofplots));
			 disp(sprintf('   figurenumber: %i',opt.figurenumber));
			 if ~isempty(opt.list),
				 disp(sprintf('   list: (%ix%i)',size(opt.list,1),size(opt.list,2)));
				 for i=1:size(opt.list,1),
					 unit=opt.list{i};
					 disp(sprintf('\n   options of plot number %i',i));
					 for j=1:size(unit.list,1)
						 if ischar(unit.list{j,2}),
							 disp(sprintf('     field: %-10s value: ''%s''',unit.list{j,1},unit.list{j,2}));
						 elseif isnumeric(unit.list{j,2}) & length(unit.list{j,2})==1,
							 disp(sprintf('     field: %-10s value: %g',unit.list{j,1},unit.list{j,2}));
						 else
							 disp(sprintf('     field: %-10s value: (%ix%i)',unit.list{j,1},size(unit.list{j,2},1),size(unit.list{j,2},2)));
						 end
					 end
				 end
			 else
				 disp(sprintf('   list: empty'));
			 end
		 end
		 %}}}
		 function opt=buildlist(opt,varargin) % {{{

			 %check length of input
			 if mod((nargin-1),2),
				 for i=1:2:(nargin-1)
					 if ~ischar(varargin{i}),
						 disp(['Last valid option: ' varargin{i-2} ]);
						 break;
					 end
				 end
				 error('Invalid parameter/value pair arguments')
			 end

			 %go through varargin and build list (like pairoptions)
			 rawoptions=pairoptions(varargin{:});
			 numoptions = (nargin-1)/2;
			 rawlist=cell(numoptions,2);
			 counter=1;
			 for i=1:numoptions,
				 optionname = varargin{2*i-1};
				 optionval  = varargin{2*i};
				 if ischar(optionname),
					 rawlist{counter,1}=optionname;
					 rawlist{counter,2}=optionval;
					 counter=counter+1;
				 else
					 %option is not a string, ignore it
					 disp(['WARNING: option number ' num2str(i) ' is not a string, it will be ignored']);
					 rawlist(counter,:)=[];
					 continue
				 end
			 end

			 %get number of data to be plotted
			 numberofplots=fieldoccurrences(rawoptions,'data');
			 opt.numberofplots=numberofplots;

			 %figure out wether alloptions flog is on
			 if strcmpi(getfieldvalue(rawoptions,'alloptions','off'),'on'),
				 allflag=1;
			 else
				 allflag=0;
			 end

			 %initialize opt.list
			 opt.list=cell(numberofplots,1);
			 for i=1:numberofplots,
				 opt.list{i}=pairoptions;
			 end

			 %process plot options
			 for i=1:size(rawlist,1),

				 %If alloptions flag has is on, apply to all plots
				 if (allflag & ~strcmpi(rawlist{i,1},'data') & ~ismember('#',rawlist{i,1})),
					 for j=1:numberofplots,
						 opt.list{j}=addfield(opt.list{j},rawlist{i,1},rawlist{i,2});
					 end

					 %option contains '#'
				 elseif ismember('#',rawlist{i,1}),

					 %get suplot(s) associated
					 string=strsplit_strict(rawlist{i,1},'#');
					 plotnums=string{end};
					 field=string{1};

					 %divide plotnums if there is a comma ','
					 plotnums=strsplit_strict(plotnums,',');

					 %loop over plotnums
					 for k=1:length(plotnums);
						 plotnum=plotnums{k};

						 %Empty
						 if isempty(plotnum),
							 continue;

							 %pound all
						 elseif strcmpi(plotnum,'all');
							 for j=1:numberofplots,
								 opt.list{j}=addfield(opt.list{j},field,rawlist{i,2});
							 end

							 %pound i-j
						 elseif ismember('-',plotnum)
							 nums=strsplit_strict(plotnum,'-');
							 if length(nums)~=2, continue; end
							 if isempty(str2num(nums{1}))|isempty(str2num(nums{2}))
								 error(['the option #i-j is not set properly for ' field]);
							 end
							 for j=str2num(nums{1}):str2num(nums{2}),
								 opt.list{j}=addfield(opt.list{j},field,rawlist{i,2});
							 end

							 %pound i
						 else
							 %assign to subplot
							 if str2num(plotnum)>numberofplots,
								 error(['opt error message: ' field ' cannot be assigned (' plotnum ' exceed maximum number of plot)']);
							 end
							 opt.list{str2num(plotnum)}=addfield(opt.list{str2num(plotnum)},field,rawlist{i,2});
						 end
					 end

					 %assign option field to corresponding subplot
				 else

					 %go through all subplot and assign to the first one free
					 j=1;
					 while (j<=numberofplots),
						 if ~exist(opt.list{j},rawlist{i,1});
							 opt.list{j}=addfield(opt.list{j},rawlist{i,1},rawlist{i,2});
							 break
						 else
							 j=j+1;
						 end
					 end
					 if j>numberofplots,
						 disp(['plot info message: too many ''' rawlist{i,1} ''' options']);
					 end
				 end
			 end

			 %check that there is no duplicates
			 for i=1:numberofplots,
				 opt.list{i}=deleteduplicates(opt.list{i},1);
			 end

			 %Get figure number (should be in options for subplot 1)
			 opt.figurenumber=getfieldvalue(opt.list{1},'figure',1);
			 removefield(opt.list{1},'figure',0);

		 end
		 %}}}
	end
end
