%MASSCONAXPBY class definition
%
%   Usage:
%      massconaxpby=massconaxpby();
%      massconaxpby=massconaxpby('name','MassCon58+35','namex','MassCon58','alpha',.5,'namey','MassCon35','beta',.5,'definitionstring','Outputdefinition1'); 
% 
%   where name is the name of the massconaxpby object, namex is the name of the first masscon, namey the name of the second masscon and alpha,beta the 
%         multiplicators. The meaning of axpby here is: resulting masscon is the linear combination (alpha *x + beta * y) 
%         of two masscons.
%
%   See also: MASSCON, REGIONALOUTPUT

classdef massconaxpby
	properties (SetAccess=public)
		%masscon axpby
		name               = '';
		definitionstring   = ''; %String that identifies this output definition uniquely, from 'Outputdefinition[1-10]'
		namex              = '';
		namey              = '';
		alpha              = NaN;
		beta               = NaN;
	end
	
	methods
		function self = extrude(self,md) % {{{
		end % }}}
		function self = massconaxpby(varargin) % {{{
			if nargin==0,
				self=setdefaultparameters(self);
			else
				%use provided options to change fields
				options=pairoptions(varargin{:});

				%get names
				self.name=getfieldvalue(options,'name','');
				self.definitionstring=getfieldvalue(options,'definitionstring');
				self.namex=getfieldvalue(options,'namex');
				self.namey=getfieldvalue(options,'namey');

				%get multiplicators: 
				self.alpha=getfieldvalue(options,'alpha');
				self.beta=getfieldvalue(options,'beta');


			end
		end % }}}
		function self = setdefaultparameters(self) % {{{
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			if ~ischar(self.name), error('masscon error message: ''name'' field should be a string!'); end
			if ~ischar(self.namex), error('masscon error message: ''namex'' field should be a string!'); end
			if ~ischar(self.namey), error('masscon error message: ''namey'' field should be a string!'); end
		
			% Create output definition cell array for check
			OutputdefinitionStringArray={};
			for i=1:100
				OutputdefinitionStringArray{i}=strcat('Outputdefinition',num2str(i));
			end

			md = checkfield(md,'fieldname','self.definitionstring','field',self.definitionstring,'values',OutputdefinitionStringArray);

			md = checkfield(md,'fieldname','self.alpha','field',self.alpha,'NaN',1,'Inf',1,'size',[1 1]);
			md = checkfield(md,'fieldname','self.betaa','field',self.beta,'NaN',1,'Inf',1,'size',[1 1]);

		end % }}}
		function md = disp(self) % {{{
		
			disp(sprintf('   Massconaxpby:\n'));

			fielddisplay(self,'name','name');
			fielddisplay(self,'definitionstring','string that identifies this output definition uniquely, from Outputdefinition[1-10]');
			fielddisplay(self,'namex','identifier for the first masscon');
			fielddisplay(self,'alpha','first masscon multiplicator');
			fielddisplay(self,'namey','identifier for the second masscon');
			fielddisplay(self,'beta','second masscon multiplicator');

		end % }}}
		function md = marshall(self,prefix,md,fid) % {{{

			WriteData(fid,prefix,'data',self.name,'name','md.massconaxpby.name','format','String');
			WriteData(fid,prefix,'data',self.definitionstring,'name','md.massconaxpby.definitionstring','format','String');
			WriteData(fid,prefix,'data',self.namex,'name','md.massconaxpby.namex','format','String');
			WriteData(fid,prefix,'data',self.namey,'name','md.massconaxpby.namey','format','String');
			WriteData(fid,prefix,'data',self.alpha,'name','md.massconaxpby.alpha','format','Double');
			WriteData(fid,prefix,'data',self.beta,'name','md.massconaxpby.beta','format','Double');

		end % }}}
	end
end
