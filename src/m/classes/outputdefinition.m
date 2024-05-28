%CONSTANTS class definition
%
%   Usage:
%      outputdefinition=outputdefinition();

classdef outputdefinition
	properties (SetAccess=public) 
		definitions                 = {};
	end
	methods
		function self = extrude(self,md) % {{{
			for i=1:length(self.definitions)
				self.definitions{i}=extrude(self.definitions{i},md);
			end
		end % }}}
		function self = outputdefinition(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = setdefaultparameters(self) % {{{

			self.definitions={};

		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			md = checkfield(md,'fieldname','outputdefinition.definitions','cell',1);

			for i=1:length(self.definitions),
				md=checkconsistency(self.definitions{i},md,solution,analyses);
			end

		end % }}}
		function disp(self) % {{{
			disp(sprintf('   outputdefinition:'));
			fielddisplay(self,'definitions','list of potential outputs that can be requested, but which need additional data to be defined');

		end % }}}
		function marshall(self,prefix,md,fid) % {{{
		data={};
		for i=1:length(self.definitions)
			self.definitions{i}.marshall(prefix,md,fid);
			classdefinition=class(self.definitions{i});
			classdefinition(1)=upper(classdefinition(1));
			data{i}=classdefinition;
		end
		data_unique=unique(data);
		WriteData(fid,prefix,'data',data_unique,'name','md.outputdefinition.list','format','StringArray');
		end % }}}
		function savemodeljs(self,fid,modelname) % {{{

			fprintf(fid,'%s.outputdefinition.definitions=[];\n',modelname);
			if ~isempty(self.definitions),
				error('outputdefinition savemodeljs error message: not supported yet!');
			end

		end % }}}
	end
end
