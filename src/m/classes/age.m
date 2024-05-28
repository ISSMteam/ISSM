%AGE class definition
%
%   Usage:
%      age=age();

classdef age
	properties (SetAccess=public)
		 spcage                 = NaN;
		 stabilization          = 0;
		 requested_outputs      = {};
	end
	methods (Static)
		function self = loadobj(self) % {{{
			% This function is directly called by matlab when a model object is
			% loaded. If the input is a struct it is an old version of this class and
			% old fields must be recovered (make sure they are in the deprecated
			% model properties)

			if verLessThan('matlab','7.9'),
				disp('Warning: your matlab version is old and there is a risk that load does not work correctly');
				disp('         if the model is not loaded correctly, rename temporarily loadobj so that matlab does not use it');

				% This is a Matlab bug: all the fields of md have their default value
				% Example of error message:
				% Warning: Error loading an object of class 'model':
				% Undefined function or method 'exist' for input arguments of type 'cell'
				%
				% This has been fixed in MATLAB 7.9 (R2009b) and later versions
			end

			if isstruct(self)
				disp('Recovering age from older version');
				self = structtoobj(age(),self);
			end
		end% }}}
	end
	methods
		function self = age(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				case 1
					inputstruct=varargin{1};
					list1 = properties('age');
					list2 = fieldnames(inputstruct);
					for i=1:length(list1)
						fieldname = list1{i};
						if ismember(fieldname,list2),
							self.(fieldname) = inputstruct.(fieldname);
						end
					end
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = extrude(self,md) % {{{
			self.spcage=project3d(md,'vector',self.spcage,'type','node');
		end % }}}
		function list = defaultoutputs(self,md) % {{{

			list = {'Age'};

		end % }}}
		function self = setdefaultparameters(self) % {{{

			%Type of stabilization to use 0:nothing 1:artificial_diffusivity
			self.stabilization=2;

			%default output
			self.requested_outputs={'default'};
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			%Early return
			if ~ismember('AgeAnalysis',analyses) |  (strcmp(solution,'TransientSolution') & md.transient.isage==0), return; end
			if dimension(md.mesh)~=3
				md = checkmessage(md,['age model only supported in 3D']);
			end

			md = checkfield(md,'fieldname','age.spcage','Inf',1,'timeseries',1);
			md = checkfield(md,'fieldname','age.stabilization','values',[0 1 2 4 5]);
			md = checkfield(md,'fieldname','age.requested_outputs','stringrow',1);
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   Masstransport solution parameters:'));
			fielddisplay(self,'spcage','Age constraint (NaN means no constraint) [yr]');
			fielddisplay(self,'stabilization','0: no stabilization, 1: artificial diffusion, 2: streamline upwinding, 4: flux corrected transport, 5: streamline upwind Petrov-Galerkin (SUPG)');
			fielddisplay(self,'requested_outputs','additional outputs requested');

		end % }}}
		function marshall(self,prefix,md,fid) % {{{

			yts=md.constants.yts;

			WriteData(fid,prefix,'object',self,'fieldname','spcage','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'fieldname','stabilization','format','Integer');

			%process requested outputs
			outputs = self.requested_outputs;
			pos  = find(ismember(outputs,'default'));
			if ~isempty(pos),
				outputs(pos) = [];                         %remove 'default' from outputs
				outputs      = [outputs defaultoutputs(self,md)]; %add defaults
			end
			WriteData(fid,prefix,'data',outputs,'name','md.age.requested_outputs','format','StringArray');
		end % }}}
		function savemodeljs(self,fid,modelname) % {{{

			writejs1Darray(fid,[modelname '.age.spcage'],self.spcage);
			writejsdouble(fid,[modelname '.age.stabilization'],self.stabilization);
			writejscellstring(fid,[modelname '.age.requested_outputs'],self.requested_outputs);

		end % }}}
	end
end
