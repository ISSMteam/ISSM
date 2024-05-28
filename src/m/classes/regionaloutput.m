%REGIONALOUTPUT class definition
%
%   Usage:
%      regionaloutput=regionaloutput();
%      regionaloutput=regionaloutput('name','Volume1','definitionstring','Outputdefinition1','outputnamestring','IceVolume',...
%                    'mask',mask);
%      regionaloutput=regionaloutput('name','Volume1','definitionstring','Outputdefinition1','outputnamestring','IceVolume',...
%                    'maskexpstring','Exp/Mask.exp','model',md)
% 
%   where mask is a vectorial field of size md.mesh.numberofvertices,1 : where vertices with values > 0 are to be included in the calculated region.
%   Alternatively, the user can pass in an Argus file and model object instead of a mask, and mask will be calculated for the user
%
%   See also: MISFIT, MASSCON, MASSCONAXPBY

classdef regionaloutput
	properties (SetAccess=public)
		%regionaloutput
		name              = '';
		definitionstring  = ''; %string that identifies this output definition uniquely, from 'Outputdefinition[1-100]'
		outputnamestring  = ''; %string that identifies the type of output you want, eg. IceVolume, TotalSmb, GroudedArea
		mask              = NaN; %mask vectorial field which identifies the region of interest (value > 0 will be included)
		maskexpstring     = '';
	end
	
	methods
		function self = extrude(self,md) % {{{
			if ~isnan(self.mask)
				self.mask=project3d(md,'vector',self.mask,'type','node');
			end
		end % }}}
		function self = regionaloutput(varargin) % {{{
			if nargin==0,
				self=setdefaultparameters(self);
			else
				%use provided options to change fields
				options=pairoptions(varargin{:});

				%get name
				self.name=getfieldvalue(options,'name','');
				if nargin==8
					self.mask=getfieldvalue(options,'mask',NaN);
					if isnan(self.mask)
						error('regionaloutput error message: ''mask'' field or ''maskexpstring'' and ''model'' fields should be defined!');
					end
				elseif nargin==10
					modelname=getfieldvalue(options,'model');
					self.maskexpstring=getfieldvalue(options,'maskexpstring');
					self=setmaskfromexp(self,modelname);
				else
					error('regionaloutput error message: ''mask'' field or ''maskexpstring'' and ''model'' fields should be defined!');
				end
		
				self.definitionstring=getfieldvalue(options,'definitionstring');
				self.outputnamestring=getfieldvalue(options,'outputnamestring');

			end
		end % }}}
		function self = setdefaultparameters(self) % {{{
		end % }}}
		function self = setmaskfromexp(self,md) % {{{

			if length(self.maskexpstring)>0
				self.mask=ContourToMesh(md.mesh.elements,md.mesh.x,md.mesh.y,self.maskexpstring,'node',1);
			end
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			if ~ischar(self.name),
				error('regionaloutput error message: ''name'' field should be a string!');
			end
			if ~ischar(self.outputnamestring),
				error('regionaloutput error message: ''outputnamestring'' field should be a string!');
			end
			
			OutputdefinitionStringArray={};
			for i=1:100
				OutputdefinitionStringArray{i}=strcat('Outputdefinition',num2str(i));
			end
			self=setmaskfromexp(self,md);

			md = checkfield(md,'fieldname','self.definitionstring','field',self.definitionstring,'values',OutputdefinitionStringArray);
			md = checkfield(md,'fieldname','self.mask','field',self.mask,'size',[md.mesh.numberofvertices 1],'NaN',1,'Inf',1);

		end % }}}
		function md = disp(self) % {{{
		
			disp(sprintf('   Regionaloutput:\n'));

			fielddisplay(self,'name','identifier for this regional response');
			fielddisplay(self,'definitionstring','string that identifies this output definition uniquely, from ''Outputdefinition[1-100]''');
			fielddisplay(self,'outputnamestring','string that identifies the type of output you want, eg. IceVolume, TotalSmb, GroudedArea');
			fielddisplay(self,'mask','mask vectorial field which identifies the region of interest (value > 0 will be included)');
			fielddisplay(self,'maskexpstring','name of Argus file that can be passed in to define the regional mask');

		end % }}}
		function md = marshall(self,prefix,md,fid) % {{{

			self=setmaskfromexp(self,md);
			WriteData(fid,prefix,'data',self.name,'name','md.regionaloutput.name','format','String');
			WriteData(fid,prefix,'data',self.definitionstring,'name','md.regionaloutput.definitionstring','format','String');
			WriteData(fid,prefix,'data',self.outputnamestring,'name','md.regionaloutput.outputnamestring','format','String');
			WriteData(fid,prefix,'data',self.mask,'name','md.regionaloutput.mask','format','DoubleMat','mattype',1);

		end % }}}
	end
end
