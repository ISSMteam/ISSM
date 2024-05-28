%MASSTRANSPORT class definition
%
%   Usage:
%      masstransport=masstransport();

classdef masstransport
	properties (SetAccess=public)
		 spcthickness           = NaN;
		 isfreesurface          = 0;
		 min_thickness          = 0;
		 hydrostatic_adjustment = 0;
		 stabilization          = 0;
		 vertex_pairing         = NaN;
		 penalty_factor         = 0;
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
				disp('Recovering masstransport from older version');
				self = structtoobj(masstransport(),self);
			end
		end% }}}
	end
	methods
		function self = masstransport(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				case 1
					inputstruct=varargin{1};
					list1 = properties('masstransport');
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
			self.spcthickness=project3d(md,'vector',self.spcthickness,'type','node');
		end % }}}
		function list = defaultoutputs(self,md) % {{{

			list = {'Thickness','Surface','Base'};

		end % }}}
		function self = setdefaultparameters(self) % {{{

			%Type of stabilization to use 0:nothing 1:artificial_diffusivity 3:Discontinuous Galerkin
			self.stabilization=1;

			%Factor applied to compute the penalties kappa=max(stiffness matrix)*10^penalty_factor
			self.penalty_factor=3;

			%Minimum ice thickness that can be used
			self.min_thickness=1;

			%Hydrostatic adjustment
			self.hydrostatic_adjustment='Absolute';

			%default output
			self.requested_outputs={'default'};
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			%Early return
			if ~ismember('MasstransportAnalysis',analyses) |  (strcmp(solution,'TransientSolution') & md.transient.ismasstransport==0), return; end

			md = checkfield(md,'fieldname','masstransport.spcthickness','Inf',1,'timeseries',1);
			md = checkfield(md,'fieldname','masstransport.isfreesurface','values',[0 1]);
			md = checkfield(md,'fieldname','masstransport.hydrostatic_adjustment','values',{'Absolute' 'Incremental'});
			md = checkfield(md,'fieldname','masstransport.stabilization','values',[0 1 2 3 4 5]);
			md = checkfield(md,'fieldname','masstransport.min_thickness','>',0);
			md = checkfield(md,'fieldname','masstransport.requested_outputs','stringrow',1);
			if ~any(isnan(md.stressbalance.vertex_pairing)),
				md = checkfield(md,'fieldname','stressbalance.vertex_pairing','>',0);
			end
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   Masstransport solution parameters:'));
			fielddisplay(self,'spcthickness','thickness constraints (NaN means no constraint) [m]');
			fielddisplay(self,'isfreesurface','do we use free surfaces (FS only) or mass conservation');
			fielddisplay(self,'min_thickness','minimum ice thickness allowed [m]');
			fielddisplay(self,'hydrostatic_adjustment','adjustment of ice shelves surface and bed elevations: ''Incremental'' or ''Absolute'' ');
			fielddisplay(self,'stabilization','0: no stabilization, 1: artificial diffusion, 2: streamline upwinding, 3: discontinuous Galerkin, 4: flux corrected transport, 5: streamline upwind Petrov-Galerkin (SUPG)');

			disp(sprintf('\n      %s','Penalty options:'));
			fielddisplay(self,'penalty_factor','offset used by penalties: penalty = Kmax*10^offset');
			fielddisplay(self,'vertex_pairing','pairs of vertices that are penalized');
			fielddisplay(self,'requested_outputs','additional outputs requested');

		end % }}}
		function marshall(self,prefix,md,fid) % {{{

			yts=md.constants.yts;

			WriteData(fid,prefix,'object',self,'fieldname','spcthickness','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'fieldname','isfreesurface','format','Boolean');
			WriteData(fid,prefix,'object',self,'fieldname','min_thickness','format','Double');
			WriteData(fid,prefix,'data',self.hydrostatic_adjustment,'format','String','name','md.masstransport.hydrostatic_adjustment');
			WriteData(fid,prefix,'object',self,'fieldname','stabilization','format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','vertex_pairing','format','DoubleMat','mattype',3);
			WriteData(fid,prefix,'object',self,'fieldname','penalty_factor','format','Double');

			%process requested outputs
			outputs = self.requested_outputs;
			pos  = find(ismember(outputs,'default'));
			if ~isempty(pos),
				outputs(pos) = [];                         %remove 'default' from outputs
				outputs      = [outputs defaultoutputs(self,md)]; %add defaults
			end
			WriteData(fid,prefix,'data',outputs,'name','md.masstransport.requested_outputs','format','StringArray');
		end % }}}
		function savemodeljs(self,fid,modelname) % {{{

			writejs1Darray(fid,[modelname '.masstransport.spcthickness'],self.spcthickness);
			writejsdouble(fid,[modelname '.masstransport.isfreesurface'],self.isfreesurface);
			writejsdouble(fid,[modelname '.masstransport.min_thickness'],self.min_thickness);
			writejsstring(fid,[modelname '.masstransport.hydrostatic_adjustment'],self.hydrostatic_adjustment);
			writejsdouble(fid,[modelname '.masstransport.stabilization'],self.stabilization);
			writejs2Darray(fid,[modelname '.masstransport.vertex_pairing'],self.vertex_pairing);
			writejsdouble(fid,[modelname '.masstransport.penalty_factor'],self.penalty_factor);
			writejscellstring(fid,[modelname '.masstransport.requested_outputs'],self.requested_outputs);

		end % }}}
	end
end
