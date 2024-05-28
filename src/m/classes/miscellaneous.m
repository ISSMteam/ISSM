%MISCELLANEOUS class definition
%
%   Usage:
%      miscellaneous=miscellaneous();

classdef miscellaneous
	properties (SetAccess=public) 
		notes = '';
		name  = '';
		dummy = struct();
	end
	methods
		function self = miscellaneous(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			md = checkfield(md,'fieldname','miscellaneous.name','empty',1);

		end % }}}
		function self = setdefaultparameters(self) % {{{

			%Add some information about the model for future reference
			username = getenv('USER');
			issmver = num2str(issmversion());
			today   = date();
			host    =  oshostname();
			self.notes = [' Model created on ' today ' by ' username ' on ' host sprintf('\n') ...
				' ISSM version: ' issmver sprintf('\n') ...
				' (path: ' pwd() ')'];

		end % }}}
		function disp(self) % {{{
			disp(sprintf('   miscellaneous parameters:'));

			fielddisplay(self,'notes','notes in a cell of strings');
			fielddisplay(self,'name','model name');
			fielddisplay(self,'dummy','empty field to store some data');

		end % }}}
		function marshall(self,prefix,md,fid) % {{{
			WriteData(fid,prefix,'object',self,'fieldname','name','format','String');
		end % }}}
		function savemodeljs(self,fid,modelname) % {{{

			writejsstring(fid,[modelname '.miscellaneous.notes'],self.notes);
			writejsstring(fid,[modelname '.miscellaneous.name'],self.name);
			if strcmpi(class(self.dummy),'double'),
				if size(self.dummy,2)==1,
					writejs1Darray(fid,[modelname '.miscellaneous.dummy'],self.dummy);
				else
					writejs2Darray(fid,[modelname '.miscellaneous.dummy'],self.dummy);
				end
			end

		end % }}}
	end
end
