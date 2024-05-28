%DEPENDENT class definition
%
%   Usage:
%      dependent=dependent();

classdef dependent
	properties (SetAccess=public) 
		name                 = '';
		fos_reverse_index    = NaN;
		exp                  = '';
		segments             = [];
		nods                 = 0;
	end
	methods
		function self = dependent(varargin) % {{{

			%use provided options to change fields
			options=pairoptions(varargin{:});

			self.name=getfieldvalue(options,'name','');
			self.exp=getfieldvalue(options,'exp','');
			self.segments=getfieldvalue(options,'segments',[]);
			self.nods=getfieldvalue(options,'nods',0);

			%if name is mass flux: 
			if strcmpi(self.name,'MassFlux'),
				%make sure that we supplied a file and that it exists! 
				if exist(self.exp)~=2,
					error('dependent checkconsistency: specified ''exp'' file does not exist!');
				end
				%process the file and retrieve segments
				mesh=getfieldvalue(options,'mesh');
				self.segments=MeshProfileIntersection(mesh.elements,mesh.x,mesh.y,self.exp);
			end
		end
		%}}}
		function self = setdefaultparameters(self) % {{{
			%do nothing
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{
			if strcmpi(self.name,'MassFlux'),
				if isempty(self.segments),
					error('dependent checkconsistency error: need segments to compute this dependent response');
				end
				if self.index<=0,
					error('dependent checkconsistency error: index for segments should be >=1');
				end
			end
			if ~isnan(self.fos_reverse_index),
				if ~strcmpi(driver,'fos_reverse'),
					error('cannot declare a dependent with a fos_reverse_index when the driver is not fos_reverse!');
				end
				if self.nods==0,
					error('dependent checkconsistency error: nods should be set to the size of the independent variable');
				end
			end

		end % }}}
		function disp(self) % {{{
			disp(sprintf('   dependent variable:'));

			fielddisplay(self,'name','variable name (must match corresponding String)');
			if ~isnan(self.fos_reverse_index),
				fielddisplay(self,'fos_reverse_index','index for fos_reverse driver of ADOLC');
			end
			if ~isempty(self.exp),
				fielddisplay(self,'exp','file needed to compute dependent variable');
				fielddisplay(self,'segments','mass flux segments');
			end

		end % }}}
	end
end
