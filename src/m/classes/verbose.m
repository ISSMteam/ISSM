%VERBOSE class definition
%
%   Available verbosity levels:
%      mprocessor  : model processing 
%      module      : modules
%      solution    : solution sequence
%      solver      : solver info (extensive)
%      convergence : convergence criteria
%      control     : control method
%      qmu         : sensitivity analysis
%      autodiff    : AD analysis
%      smb         : smb analysis
%
%   Usage:
%      verbose=verbose();
%      verbose=verbose(3);
%      verbose=verbose('all');
%      verbose=verbose('001100');
%      verbose=verbose('module',true,'solver',false);

%WARNING: some parts of this file are Synchronized with src/c/shared/Numerics/Verbosity.h
%         Do not modify these sections. See src/c/shared/Numerics/README for more info

classdef verbose
	properties (SetAccess=public) 
		% {{{
		%BEGINFIELDS
		mprocessor=false;
		module=false;
		solution=false;
		solver=false;
		convergence=false;
		control=false;
		qmu=false;
		autodiff=false;
		smb=false;
		%ENDFIELDS
		% }}}
	end
	%}}}
	methods
		function verbose=verbose(varargin) % {{{

			switch(nargin),
				case 0,
					verbose.solution=true;
					verbose.qmu=true;
					verbose.control=true;
				case 1,
					binary=varargin{1};
					if     ischar(binary),
						if strcmpi(binary,'all'),
							binary=2^11-1; %all ones
							verbose=BinaryToVerbose(verbose,binary);
							verbose.solver=false; %Do not use by default
						else
							binary=bin2dec(binary);
							verbose=BinaryToVerbose(verbose,binary);
						end
					elseif isnumeric(binary),
						verbose=BinaryToVerbose(verbose,binary);
					end 
				otherwise,
					%Use options to initialize object
					verbose=AssignObjectFields(pairoptions(varargin{:}),verbose);

					%Cast to logicals
					listproperties=properties('verbose');
					for i=1:numel(listproperties),
						fieldname=listproperties{i};
						fieldvalue=verbose.(fieldname);
						if (islogical(fieldvalue) | isnumeric(fieldvalue)) & numel(fieldvalue)==1,
							verbose.(fieldname)=logical(fieldvalue);
						else
							error('verbose supported field values are logicals only (true or false)');
						end
					end
			end
		end
		%}}}
		function binary=VerboseToBinary(verbose) % {{{

		%BEGINVERB2BIN
		binary=0;
		if (verbose.mprocessor), binary=bitor(binary,1); end
		if (verbose.module), binary=bitor(binary,2); end
		if (verbose.solution), binary=bitor(binary,4); end
		if (verbose.solver), binary=bitor(binary,8); end
		if (verbose.convergence), binary=bitor(binary,16); end
		if (verbose.control), binary=bitor(binary,32); end
		if (verbose.qmu), binary=bitor(binary,64); end
		if (verbose.autodiff), binary=bitor(binary,128); end
		if (verbose.smb), binary=bitor(binary,256); end
		%ENDVERB2BIN

		end
		%}}}
		function verbose=BinaryToVerbose(verbose,binary) % {{{

		%BEGINBIN2VERB
		if bitand(binary,1), verbose.mprocessor=true; else verbose.mprocessor=false; end
		if bitand(binary,2), verbose.module=true; else verbose.module=false; end
		if bitand(binary,4), verbose.solution=true; else verbose.solution=false; end
		if bitand(binary,8), verbose.solver=true; else verbose.solver=false; end
		if bitand(binary,16), verbose.convergence=true; else verbose.convergence=false; end
		if bitand(binary,32), verbose.control=true; else verbose.control=false; end
		if bitand(binary,64), verbose.qmu=true; else verbose.qmu=false; end
		if bitand(binary,128), verbose.autodiff=true; else verbose.autodiff=false; end
		if bitand(binary,256), verbose.smb=true; else verbose.smb=false; end
		%ENDBIN2VERB

		end
		%}}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			if md.inversion.iscontrol
				temp = verbose('control',1);
				if(VerboseToBinary(self) ~= VerboseToBinary(temp))
					disp('INFO: the outlog will look better if only md.verbose.control is turned on');
				end
			end

		end % }}}
		function disp(verbose) % {{{

		%BEGINDISP
		disp(sprintf('class ''%s''  = ',class(verbose)));
		disp(sprintf('   %15s : %s','mprocessor',mat2str(verbose.mprocessor)));
		disp(sprintf('   %15s : %s','module',mat2str(verbose.module)));
		disp(sprintf('   %15s : %s','solution',mat2str(verbose.solution)));
		disp(sprintf('   %15s : %s','solver',mat2str(verbose.solver)));
		disp(sprintf('   %15s : %s','convergence',mat2str(verbose.convergence)));
		disp(sprintf('   %15s : %s','control',mat2str(verbose.control)));
		disp(sprintf('   %15s : %s','qmu',mat2str(verbose.qmu)));
		disp(sprintf('   %15s : %s','autodiff',mat2str(verbose.autodiff)));
		disp(sprintf('   %15s : %s','smb',mat2str(verbose.smb)));
		%ENDDISP

		end
		%}}}
		function marshall(self,prefix,md,fid) % {{{
			WriteData(fid,prefix,'data',VerboseToBinary(self),'name','md.verbose','format','Integer');
		end % }}}
		function savemodeljs(self,fid,modelname) % {{{
		
			writejsdouble(fid,[modelname '.verbose.mprocessor'],self.mprocessor);
			writejsdouble(fid,[modelname '.verbose.module'],self.module);
			writejsdouble(fid,[modelname '.verbose.solution'],self.solution);
			writejsdouble(fid,[modelname '.verbose.solver'],self.solver);
			writejsdouble(fid,[modelname '.verbose.convergence'],self.convergence);
			writejsdouble(fid,[modelname '.verbose.control'],self.control);
			writejsdouble(fid,[modelname '.verbose.qmu'],self.qmu);
			writejsdouble(fid,[modelname '.verbose.autodiff'],self.autodiff);
			writejsdouble(fid,[modelname '.verbose.smb'],self.smb);

		end % }}}
	end
end
