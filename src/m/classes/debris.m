%DEBRIS class definition
%
%   Usage:
%      debris=debris();

classdef debris
	properties (SetAccess=public)
		spcthickness             = NaN;
		min_thickness            = 0;
		stabilization            = 0;
		packingfraction          = 0;
		removalmodel             = 0;
		displacementmodel        = 0;
		max_displacementvelocity = 0;
		removal_slope_threshold  = 0;
		removal_stress_threshold = 0;
		vertex_pairing           = NaN;
		requested_outputs        = {};
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
				disp('Recovering debris from older version');
				self = structtoobj(debris(),self);
			end
		end% }}}
	end
	methods
		function self = debris(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				case 1
					inputstruct=varargin{1};
					list1 = properties('debris');
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

			list = {'DebrisThickness','DebrisMaskNodeActivation','VxDebris','VyDebris'};

		end % }}}
		function self = setdefaultparameters(self) % {{{

			%Type of stabilization to use 0:nothing 1:artificial_diffusivity 3:Discontinuous Galerkin
			self.stabilization=2;

			%Minimum debris thickness that can be used
			self.min_thickness=0;

			%Fraction of debris covered in the ice
			self.packingfraction=0.01;

			%Type of frontal debris removal
			self.removalmodel=0;

			%Type of debris displacement
			self.displacementmodel=0;

			%Slope threshold for removalmodel (1)
			self.removal_slope_threshold=0;

			%Stress threshold for removalmodel (2)
			self.removal_stress_threshold=0;

			%Max velocity for displacementmodel (1)
			self.max_displacementvelocity=0;

			%default output
			self.requested_outputs={'default'};
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			%Early return
			if ~ismember('MasstransportAnalysis',analyses) | (strcmp(solution,'TransientSolution') & md.transient.isdebris==0), return; end

			md = checkfield(md,'fieldname','debris.spcthickness');
			md = checkfield(md,'fieldname','debris.stabilization','values',[0 1 2 3 4 5]);
			md = checkfield(md,'fieldname','debris.min_thickness','>=',0);
			md = checkfield(md,'fieldname','debris.packingfraction','>=',0);
			md = checkfield(md,'fieldname','debris.removalmodel','values',[0 1 2]);
			md = checkfield(md,'fieldname','debris.displacementmodel','values',[0 1 2]);
			md = checkfield(md,'fieldname','debris.max_displacementvelocity','>=',0);
			md = checkfield(md,'fieldname','debris.removal_slope_threshold','>=',0);
			md = checkfield(md,'fieldname','debris.removal_stress_threshold','>=',0);
			md = checkfield(md,'fieldname','debris.requested_outputs','stringrow',1);
			if ~any(isnan(md.stressbalance.vertex_pairing)),
				md = checkfield(md,'fieldname','stressbalance.vertex_pairing','>',0);
			end
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   debris solution parameters:'));
			fielddisplay(self,'spcthickness','debris thickness constraints (NaN means no constraint) [m]');
			fielddisplay(self,'min_thickness','minimum debris thickness allowed [m]');
			fielddisplay(self,'packingfraction','fraction of debris covered in the ice');
			fielddisplay(self,'stabilization','0: no stabilization, 1: artificial diffusion, 2: streamline upwinding, 3: streamline upwind Petrov-Galerkin (SUPG)');
			fielddisplay(self,'removalmodel','frontal removal of debris. 0: no removal, 1: Slope-triggered debris removal, 2: driving-stress triggered debris removal');
			fielddisplay(self,'displacementmodel','debris displacement. 0: no displacement, 1: additional debris velocity above the critical slope/stress threshold');
			fielddisplay(self,'max_displacementvelocity','maximum velocity of debris transport (v_ice + v_displacement) (m/a)');
			fielddisplay(self,'removal_slope_threshold','critical slope (degrees) for removalmodel (1)');
			fielddisplay(self,'removal_stress_threshold','critical stress (Pa) for removalmodel (2)');

			disp(sprintf('\n      %s','Penalty options:'));
			fielddisplay(self,'vertex_pairing','pairs of vertices that are penalized');
			fielddisplay(self,'requested_outputs','additional outputs requested');

		end % }}}
		function marshall(self,prefix,md,fid) % {{{
			WriteData(fid,prefix,'object',self,'fieldname','spcthickness','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'fieldname','min_thickness','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','stabilization','format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','removalmodel','format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','displacementmodel','format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','max_displacementvelocity','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','removal_slope_threshold','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','removal_stress_threshold','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','packingfraction','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','vertex_pairing','format','DoubleMat','mattype',3);

			%process requested outputs
			outputs = self.requested_outputs;
			pos  = find(ismember(outputs,'default'));
			if ~isempty(pos),
				outputs(pos) = [];                         %remove 'default' from outputs
				outputs      = [outputs defaultoutputs(self,md)]; %add defaults
			end
			WriteData(fid,prefix,'data',outputs,'name','md.debris.requested_outputs','format','StringArray');
		end % }}}
		function savemodeljs(self,fid,modelname) % {{{

			writejs1Darray(fid,[modelname '.debris.spcthickness'],self.spcthickness);
			writejsdouble(fid,[modelname '.debris.min_thickness'],self.min_thickness);
			writejsdouble(fid,[modelname '.debris.stabilization'],self.stabilization);
			writejsdouble(fid,[modelname '.debris.removalmodel'],self.removalmodel);
			writejsdouble(fid,[modelname '.debris.displacementmodel'],self.displacementmodel);
			writejsdouble(fid,[modelname '.debris.max_displacementvelocity'],self.displacementmodel);
			writejsdouble(fid,[modelname '.debris.removal_slope_threshold'],self.removal_slope_threshold);
			writejsdouble(fid,[modelname '.debris.removal_stress_threshold'],self.removal_stress_threshold);
			writejsdouble(fid,[modelname '.debris.packingfraction'],self.packingfraction);
			writejs2Darray(fid,[modelname '.debris.vertex_pairing'],self.vertex_pairing);
			writejsdouble(fid,[modelname '.debris.penalty_factor'],self.penalty_factor);
			writejscellstring(fid,[modelname '.debris.requested_outputs'],self.requested_outputs);

		end % }}}
	end
end
