
%Hydrologydc class definition
%
%   Usage:
%      hydrologydc=hydrologydc();

classdef hydrologydc
	properties (SetAccess=public)
		water_compressibility    = 0;
		isefficientlayer         = 0;
		penalty_factor           = 0;
		penalty_lock             = 0;
		rel_tol                  = 0;
		max_iter                 = 0;
		steps_per_step           = 0;
		step_adapt               = 0;
		averaging                = 0;
		sedimentlimit_flag       = 0;
		sedimentlimit            = 0;
		transfer_flag            = 0;
		unconfined_flag          = 0;
		leakage_factor           = 0;
		basal_moulin_input       = NaN;
		requested_outputs        = {};

		spcsediment_head         = NaN;
		mask_thawed_node         = NaN;
		sediment_transmitivity   = NaN;
		sediment_compressibility = 0;
		sediment_porosity        = 0;
		sediment_thickness       = 0;


		spcepl_head              = NaN;
		mask_eplactive_node      = NaN;
		epl_compressibility      = 0;
		epl_porosity             = 0;
		epl_initial_thickness    = 0;
		epl_colapse_thickness    = 0;
		epl_thick_comp           = 0;
		epl_max_thickness        = 0;
		epl_conductivity         = 0;
		eplflip_lock             = 0;
	end
	methods
		function self = extrude(self,md)    % {{{
			self.spcsediment_head=project3d(md,'vector',self.spcsediment_head,'type','node','layer',1);
			self.sediment_transmitivity=project3d(md,'vector',self.sediment_transmitivity,'type','node','layer',1);
			self.basal_moulin_input=project3d(md,'vector',self.basal_moulin_input,'type','node','layer',1);
			self.mask_thawed_node=project3d(md,'vector',self.mask_thawed_node,'type','node','layer',1);
			if(self.isefficientlayer==1);
				self.spcepl_head=project3d(md,'vector',self.spcepl_head,'type','node','layer',1);
				self.mask_eplactive_node=project3d(md,'vector',self.mask_eplactive_node,'type','node','layer',1);
			end
		end    % }}}
		function self = hydrologydc(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				otherwise
					error('constructor not supported');
			end
		end% }}}
		function list = defaultoutputs(self,md) % {{{
			list = {'SedimentHead','SedimentHeadResidual','EffectivePressure'};
			if self.isefficientlayer,
				list=[list,{'EplHead','HydrologydcMaskEplactiveNode','HydrologydcMaskEplactiveElt','EplHeadSlopeX','EplHeadSlopeY','HydrologydcEplThickness'}];
			end
			if self.steps_per_step>1 | self.step_adapt,
				list = [list,'EffectivePressureSubstep','SedimentHeadSubstep'];
				if self.isefficientlayer,
					list = [list,'EplHeadSubstep','HydrologydcEplThicknessSubstep'];
				end
			end
		end % }}}
		function self = initialize(self,md) % {{{
			self.epl_colapse_thickness = self.sediment_transmitivity/self.epl_conductivity;
			if isnan(self.basal_moulin_input),
				self.basal_moulin_input=zeros(md.mesh.numberofvertices,1);
				disp('      no hydrology.basal_moulin_input specified: values set as zero');
			end

		end % }}}
		function self = setdefaultparameters(self)% {{{
			%Parameters from de Fleurian 2014
			self.water_compressibility    = 5.04e-10;
			self.isefficientlayer         = 1;
			self.penalty_factor           = 3;
			self.penalty_lock             = 0;
			self.rel_tol                  = 1.0e-06;
			self.max_iter                 = 100;
			self.steps_per_step           = 1;
			self.step_adapt               = 0;
			self.averaging                = 0;
			self.sedimentlimit_flag       = 0;
			self.sedimentlimit            = 0;
			self.transfer_flag            = 1;
			self.unconfined_flag          = 0;
			self.leakage_factor           = 1.0e-10;
			self.requested_outputs        = {'default'};


			self.sediment_compressibility = 1.0e-08;
			self.sediment_porosity        = 0.4;
			self.sediment_thickness       = 20.0;
			self.sediment_transmitivity   = 8.0e-04;

			self.epl_compressibility      = 1.0e-08;
			self.epl_conductivity         = 8.0e-02;
			self.epl_porosity             = 0.4;
			self.epl_initial_thickness    = 1.0;
			self.epl_colapse_thickness    = self.sediment_transmitivity/self.epl_conductivity;
			self.epl_thick_comp           = 1;
			self.epl_max_thickness        = 5.0;
			self.eplflip_lock             = 0;
		end  % }}}
		function md = checkconsistency(self,md,solution,analyses)% {{{
		%Early return
			if ~ismember('HydrologyDCInefficientAnalysis',analyses) & ~ismember('HydrologyDCEfficientAnalysis',analyses),
				return;
			end

			md = checkfield(md,'fieldname','hydrology.water_compressibility','>',0,'numel',1);
			md = checkfield(md,'fieldname','hydrology.isefficientlayer','numel',[1],'values',[0 1]);
			md = checkfield(md,'fieldname','hydrology.penalty_factor','>',0,'numel',1);
			md = checkfield(md,'fieldname','hydrology.penalty_lock','>=',0,'numel',1);
			md = checkfield(md,'fieldname','hydrology.rel_tol','>',0,'numel',1);
			md = checkfield(md,'fieldname','hydrology.max_iter','>',0,'numel',1);
			md = checkfield(md,'fieldname','hydrology.steps_per_step','>',0,'numel',1);
			md = checkfield(md,'fieldname','hydrology.step_adapt','numel',1,'values',[0 1]);
			md = checkfield(md,'fieldname','hydrology.averaging','numel',[1],'values',[0 1 2]);
			md = checkfield(md,'fieldname','hydrology.sedimentlimit_flag','numel',[1],'values',[0 1 2 3]);
			md = checkfield(md,'fieldname','hydrology.transfer_flag','numel',[1],'values',[0 1]);
			md = checkfield(md,'fieldname','hydrology.unconfined_flag','numel',[1],'values',[0 1]);
			md = checkfield(md,'fieldname','hydrology.requested_outputs','stringrow',1);
			if self.sedimentlimit_flag==1,
				md = checkfield(md,'fieldname','hydrology.sedimentlimit','>',0,'numel',1);
			end
			if self.transfer_flag==1,
				md = checkfield(md,'fieldname','hydrology.leakage_factor','>',0,'numel',1);
			end
			md = checkfield(md,'fieldname','hydrology.basal_moulin_input','NaN',1,'Inf',1,'timeseries',1);

			md = checkfield(md,'fieldname','hydrology.spcsediment_head','Inf',1,'timeseries',1);
			md = checkfield(md,'fieldname','hydrology.sediment_compressibility','>',0,'numel',1);
			md = checkfield(md,'fieldname','hydrology.sediment_porosity','>',0,'numel',1);
			md = checkfield(md,'fieldname','hydrology.sediment_thickness','>',0,'numel',1);
			md = checkfield(md,'fieldname','hydrology.sediment_transmitivity','>=',0,'size',[md.mesh.numberofvertices 1]);
			md = checkfield(md,'fieldname','hydrology.mask_thawed_node','size',[md.mesh.numberofvertices 1],'values',[0 1]);

			if self.isefficientlayer==1,
				md = checkfield(md,'fieldname','hydrology.spcepl_head','Inf',1,'timeseries',1);
				md = checkfield(md,'fieldname','hydrology.mask_eplactive_node','size',[md.mesh.numberofvertices 1],'values',[0 1]);
				md = checkfield(md,'fieldname','hydrology.epl_compressibility','>',0,'numel',1);
				md = checkfield(md,'fieldname','hydrology.epl_porosity','>',0,'numel',1);
				md = checkfield(md,'fieldname','hydrology.epl_initial_thickness','>',0,'numel',1);
				md = checkfield(md,'fieldname','hydrology.epl_colapse_thickness','>',0,'numel',1);
				md = checkfield(md,'fieldname','hydrology.epl_thick_comp','numel',[1],'values',[0 1]);
				md = checkfield(md,'fieldname','hydrology.epl_max_thickness','>',0,'numel',1);
				md = checkfield(md,'fieldname','hydrology.epl_conductivity','>',0,'numel',1);
				md = checkfield(md,'fieldname','hydrology.eplflip_lock','>=',0,'numel',1);
				if (self.epl_colapse_thickness>self.epl_initial_thickness),
					md = checkmessage(md,'Colapsing thickness for EPL larger than initial thickness');
				end
			end
		end% }}}
		function disp(self)% {{{
			disp(sprintf('   hydrology Dual Porous Continuum Equivalent parameters:'));
			disp(sprintf('   - general parameters'));
			fielddisplay(self,'water_compressibility','compressibility of water [Pa^-1]');
			fielddisplay(self,'isefficientlayer','do we use an efficient drainage system [1: true; 0: false]');
			fielddisplay(self,'penalty_factor','exponent of the value used in the penalisation method [dimensionless]');
			fielddisplay(self,'penalty_lock','stabilize unstable constraints that keep zigzagging after n iteration (default is 0, no stabilization)');
			fielddisplay(self,'rel_tol','tolerance of the nonlinear iteration for the transfer between layers [dimensionless]');
			fielddisplay(self,'max_iter','maximum number of nonlinear iteration');
			fielddisplay(self,'steps_per_step','number of hydrology steps per timestep');
			fielddisplay(self,'step_adapt', 'adaptative sub stepping  [1: true 0: false] default is 0');
			fielddisplay(self, 'averaging', 'averaging methods from short to long steps');
			disp(sprintf('%55s  0: Arithmetic (default)'));
			disp(sprintf('%55s  0: Geometric'));
			disp(sprintf('%55s  0: Harmonic'));
			fielddisplay(self,'sedimentlimit_flag','what kind of upper limit is applied for the inefficient layer');
			disp(sprintf('%55s  0: no limit',' '));
			disp(sprintf('%55s  1: user defined: %s',' ','sedimentlimit'));
			disp(sprintf('%55s  2: hydrostatic pressure',' '));
			disp(sprintf('%55s  3: normal stress',' '));
			if self.sedimentlimit_flag==1,
				fielddisplay(self,'sedimentlimit','user defined upper limit for the inefficient layer [m]');
			end
			fielddisplay(self,'transfer_flag','what kind of transfer method is applied between the layers');
			disp(sprintf('%55s  0: no transfer',' '));
			disp(sprintf('%55s  1: constant leakage factor: %s',' ','leakage_factor'));
			if self.transfer_flag==1,
				fielddisplay(self,'leakage_factor','user defined leakage factor [m]');
			end
			fielddisplay(self,'unconfined_flag','Do you want unconfined scheme to be used (transitory)');
			disp(sprintf('%55s  0: confined only',' '));
			disp(sprintf('%55s  1: confined unconfined'));
      fielddisplay(self,'requested_outputs','additional outputs requested');
			fielddisplay(self,'basal_moulin_input','water flux at a given point [m3 s-1]');
			disp(sprintf('   - for the sediment layer'));
			fielddisplay(self,'spcsediment_head','sediment water head constraints (NaN means no constraint) [m above MSL]');
			fielddisplay(self,'sediment_compressibility','sediment compressibility [Pa^-1]');
			fielddisplay(self,'sediment_porosity','sediment [dimensionless]');
			fielddisplay(self,'sediment_thickness','sediment thickness [m]');
			fielddisplay(self,'sediment_transmitivity','sediment transmitivity [m^2/s]');
      fielddisplay(self,'mask_thawed_node','deactivate (0) hydrology on frozen nodes');

			if self.isefficientlayer==1,
				disp(sprintf('   - for the epl layer'));
				fielddisplay(self,'spcepl_head','epl water head constraints (NaN means no constraint) [m above MSL]');
				fielddisplay(self,'mask_eplactive_node','active (1) or not (0) EPL');
				fielddisplay(self,'epl_compressibility','epl compressibility [Pa^-1]');
				fielddisplay(self,'epl_porosity','epl [dimensionless]');
				fielddisplay(self,'epl_initial_thickness','epl initial thickness [m]');
				fielddisplay(self,'epl_colapse_thickness','epl colapsing thickness [m]');
				fielddisplay(self,'epl_thick_comp','epl thickness computation flag');
				fielddisplay(self,'epl_max_thickness','epl maximal thickness [m]');
				fielddisplay(self,'epl_conductivity','epl conductivity [m^2/s]');
				fielddisplay(self,'eplflip_lock','lock the epl activation to avoid fli-floping (default is 0, no stabilization)');
			end
		end % }}}
		function marshall(self,prefix,md,fid)% {{{
			WriteData(fid,prefix,'name','md.hydrology.model','data',1,'format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','water_compressibility','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','isefficientlayer','format','Boolean');
			WriteData(fid,prefix,'object',self,'fieldname','penalty_factor','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','penalty_lock','format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','rel_tol','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','max_iter','format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','steps_per_step','format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','step_adapt','format','Boolean');
			WriteData(fid,prefix,'object',self,'fieldname','averaging','format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','sedimentlimit_flag','format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','transfer_flag','format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','unconfined_flag','format','Integer');
			if self.sedimentlimit_flag==1,
				WriteData(fid,prefix,'object',self,'fieldname','sedimentlimit','format','Double');
			end
			if self.transfer_flag==1,
				WriteData(fid,prefix,'object',self,'fieldname','leakage_factor','format','Double');
			end
			WriteData(fid,prefix,'object',self,'fieldname','basal_moulin_input','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts)

			WriteData(fid,prefix,'object',self,'fieldname','spcsediment_head','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'fieldname','sediment_compressibility','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','sediment_porosity','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','sediment_thickness','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','sediment_transmitivity','format','DoubleMat','mattype',1');
			WriteData(fid,prefix,'object',self,'fieldname','mask_thawed_node','format','DoubleMat','mattype',1);
			if self.isefficientlayer==1,
				WriteData(fid,prefix,'object',self,'fieldname','spcepl_head','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
				WriteData(fid,prefix,'object',self,'fieldname','mask_eplactive_node','format','DoubleMat','mattype',1);
				WriteData(fid,prefix,'object',self,'fieldname','epl_compressibility','format','Double');
				WriteData(fid,prefix,'object',self,'fieldname','epl_porosity','format','Double');
				WriteData(fid,prefix,'object',self,'fieldname','epl_initial_thickness','format','Double');
				WriteData(fid,prefix,'object',self,'fieldname','epl_colapse_thickness','format','Double');
				WriteData(fid,prefix,'object',self,'fieldname','epl_thick_comp','format','Integer');
				WriteData(fid,prefix,'object',self,'fieldname','epl_max_thickness','format','Double');
				WriteData(fid,prefix,'object',self,'fieldname','epl_conductivity','format','Double');
				WriteData(fid,prefix,'object',self,'fieldname','eplflip_lock','format','Integer');
			end
			outputs = self.requested_outputs;
			pos  = find(ismember(outputs,'default'));
			if ~isempty(pos),
				outputs(pos) = [];  %remove 'default' from outputs
				outputs      = [outputs defaultoutputs(self,md)]; %add defaults
			end
			WriteData(fid,prefix,'data',outputs,'name','md.hydrology.requested_outputs','format','StringArray');
		end% }}}
	end
end
