%MASSFLUXATGATE class definition
%
%   Usage:
%      massfluxatgate=massfluxatgate();
%      massfluxatgate=massfluxatgate('name','GateName','definitionname','Outputdefinition1','profilename','PathToExpFile');

classdef massfluxatgate
	properties (SetAccess=public)  
		%massfluxatgate 
		name            = '';
		definitionstring  = '';
		profilename     = ''; 
	end
	properties (SetAccess=private)  
		segments        = NaN;
	end
	methods
		function self = extrude(self,md) % {{{
		end % }}}
		function self = massfluxatgate(varargin) % {{{
			if nargin==0,
				self=setdefaultparameters(self);
			else
				%use provided options to change fields
				options=pairoptions(varargin{:});

				%get name
				self.name=getfieldvalue(options,'name','');
				self.definitionstring=getfieldvalue(options,'definitionstring');
				self.profilename=getfieldvalue(options,'profilename');
			end
		end % }}}
		function self = setdefaultparameters(self) % {{{
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{
			
			if ~ischar(self.name),
				error('massfluxatgate error message: ''name'' field should be a string!');
			end
			if ~ischar(self.profilename),
				error('massfluxatgate error message: ''profilename'' field should be a string!');
			end
		
			OutputdefinitionStringArray={};
			for i=1:100
				OutputdefinitionStringArray{i}=strcat('Outputdefinition',num2str(i));
			end
			md = checkfield(md,'field',self.definitionstring,'values',OutputdefinitionStringArray);

			%check the profilename points to a file!: 
			if exist(self.profilename,'file')~=2,
				error('massfluxatgate error message: file name for profile corresponding to gate does not point to a legitimate file on disk!');
			end
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   Massfluxatgate:\n'));

			fielddisplay(self,'name','identifier for this massfluxatgate response');
			fielddisplay(self,'profilename','name of file (shapefile or argus file) defining a profile (or gate)');
			fielddisplay(self,'definitionstring','string that identifies this output definition uniquely, from ''Outputdefinition[1-100]''');
			
		end % }}}
		function marshall(self,prefix,md,fid) % {{{

		%before marshalling, we need to create the segments out of the profilename: 
		self.segments=MeshProfileIntersection(md.mesh.elements,md.mesh.x,md.mesh.y,self.profilename);

		%ok, marshall name and segments: 
		WriteData(fid,prefix,'data',self.name,'name','md.massfluxatgate.name','format','String');
		WriteData(fid,prefix,'data',self.definitionstring,'name','md.massfluxatgate.definitionstring','format','String');
		WriteData(fid,prefix,'data',self.segments,'name','md.massfluxatgate.segments','format','DoubleMat','mattype',1);

		end % }}}
	end
end
