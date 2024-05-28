%SOLIDEARTHSETTINGS class definition
%
%   Usage:
%      solidearthsettings=solidearthsettings();

classdef solidearthsettings
	properties (SetAccess=public)
		reltol                 = 0;
		abstol                 = 0;
		maxiter                = 0;
		selfattraction         = 1;
		elastic                = 1;
		viscous                = 1;
		rotation               = 1;
		grdocean               = 1;
		ocean_area_scaling     = 0;
		runfrequency           = 1; %how many time steps we skip before we run grd_core
		sealevelloading        = 1; %will sea-level loads be computed? 
		isgrd                  = 0; %will GRD patterns be computed? 
		compute_bp_grd         = 0; %will GRD patterns for bottom pressures be computed? 
		degacc                 = 0; %degree increment for resolution of Green tables.
		timeacc                = 1; %time step accuracy required to compute Green tables
		horiz                  = 0; %compute horizontal deformation
		grdmodel               = 1; %grd model (0 by default, 1 for (visco-)elastic, 2 for Ivins)
		cross_section_shape    = 0; %cross section only used when grd model is Ivins
	end
	methods (Static)
		function self = loadobj(self) % {{{
			% This function is directly called by matlab when a model object is
			% loaded. If the input is a struct it is an old version of this class and
			% old fields must be recovered (make sure they are in the deprecated
			% model properties)

			if isstruct(self)
				% 2021, Jun 4
				if isfield(self,'rigid')
					self.selfattraction = self.rigid;
				end
				if isfield(self,'computesealevelchange')
					self.sealevelloading = self.computesealevelchange;
				end
				self = structtoobj(solidearthsettings(),self);

			end
		end% }}}
	end
	methods

		function self = solidearthsettings(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = setdefaultparameters(self) % {{{

			%Convergence criterion: absolute, relative and residual
			self.reltol=0.01; % 1 percent
			self.abstol=NaN;  % default

			%maximum of non-linear iterations.
			self.maxiter=5;

			%computational flags:
			self.selfattraction=1;
			self.elastic=1;
			self.viscous=1;
			self.rotation=1;
			self.grdocean=1;
			self.ocean_area_scaling=0;
			self.compute_bp_grd=0;
			self.isgrd=0;
			self.sealevelloading=1;

			%numerical discretization accuracy
			self.degacc=.01;
			self.timeacc=1; 

			%how many time steps we skip before we run solidearthsettings solver during transient
			self.runfrequency=1;

			%horizontal displacement? (not on by default)
			self.horiz=0;

			%cross section for Ivins model
			self.cross_section_shape=1; %square as default (see iedge in GiaDeflectionCorex)

			%grd model by default
			self.grdmodel=1;

		end % }}}
		function disp(self) % {{{
			disp(sprintf('   solidearth settings:'));
			disp(sprintf('      core:'));
			fielddisplay(self,'isgrd','compute GRD patterns (default: 1)');
			fielddisplay(self,'grdmodel','type of deformation model, 0 for no GRD, 1 for spherical GRD model (SESAW model), 2 for half-space planar GRD (visco-elastic model from Ivins)');
			fielddisplay(self,'runfrequency','How many time steps we let masstransport core accumulate changes before each run of the sealevelchange core (default: 1, i.e run slc every time step)');
			disp(sprintf('      computational flags:'));
			fielddisplay(self,'selfattraction','enables surface mass load to perturb the gravity field');
			fielddisplay(self,'elastic','enables elastic deformation from surface loading');
			fielddisplay(self,'viscous','enables viscous deformation from surface loading');
			fielddisplay(self,'rotation','enables polar motion to feedback on the GRD fields');
			fielddisplay(self,'compute_bp_grd','compute GRD patterns for bottom pressure loads (default: 1)');
			fielddisplay(self,'cross_section_shape','1: square-edged (default). 2: elliptical. See iedge in GiaDeflectionCore. Used only for grdmodel=2 only');
			disp(sprintf('      resolution:'));
			fielddisplay(self,'degacc','spatial accuracy (default: .01 deg) for numerical discretization of the Green''s functions');
			fielddisplay(self,'timeacc','time accuracy (default: 1 yr) for numerical discretization of the Green''s functions');
			disp(sprintf('      sea-level equation:'));
			fielddisplay(self,'grdocean','does this planet have an ocean, if set to 1: global water mass is conserved in GRD module (default: 1)'); 
			fielddisplay(self,'sealevelloading','enables surface loading from sea-level change (default: 1)');
			fielddisplay(self,'maxiter','maximum number of nonlinear iterations');
			fielddisplay(self,'reltol','sea level change relative convergence criterion (default, NaN: not applied)');
			fielddisplay(self,'abstol','sea level change absolute convergence criterion(default, NaN: not applied)');
			fielddisplay(self,'ocean_area_scaling','correction for model representation of ocean area (default: No correction)'); 

		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			if ~ismember('SealevelchangeAnalysis',analyses) | (strcmp(solution,'TransientSolution') & md.transient.isslc==0), 
				return; 
			end
			md = checkfield(md,'fieldname','solidearth.settings.reltol','size',[1 1]);
			md = checkfield(md,'fieldname','solidearth.settings.abstol','size',[1 1]);
			md = checkfield(md,'fieldname','solidearth.settings.maxiter','size',[1 1],'>=',1);
			md = checkfield(md,'fieldname','solidearth.settings.runfrequency','size',[1 1],'>=',1);
			md = checkfield(md,'fieldname','solidearth.settings.degacc','size',[1 1],'>=',1e-10);
			md = checkfield(md,'fieldname','solidearth.settings.timeacc','size',[1 1],'>',0);
			md = checkfield(md,'fieldname','solidearth.settings.horiz','NaN',1,'Inf',1,'values',[0 1]);
			md = checkfield(md,'fieldname','solidearth.settings.grdmodel','>=',0,'<=',2);
			md = checkfield(md,'fieldname','solidearth.settings.cross_section_shape','numel',[1],'values',[1,2]);

			if self.elastic==1 & self.selfattraction==0,
				error('solidearthsettings checkconsistency error message: need selfattraction on if elastic flag is set');
			end
			if self.viscous==1 & self.elastic==0,
				error('solidearthsettings checkconsistency error message: need elastic on if viscous flag is set');
			end
			if self.rotation==1 & self.elastic==0,
				error('solidearthsettings checkconsistency error message: need elastic on if rotation flag is set');
			end

			%a GRD computation has been requested, make some checks on the nature of the meshes provided. 
			if self.isgrd,
				if strcmpi(class(md.mesh),'mesh3dsurface'),
					if self.grdmodel==2,
						error('model requires a 2D mesh to run gia Ivins computations (change mesh from mesh3dsurface to mesh2d)');
					end
				else
					if self.grdmodel==1,
						error('model requires a 3D surface mesh to run GRD computations (change mesh from mesh2d to mesh3dsurface)');
					end
				end
				if self.sealevelloading==1 & self.grdocean==0
					error('solidearthsettings checkconsistency error message: need grdocean on if sealevelloading flag is set');
				end
			end

			if self.compute_bp_grd==1 & md.solidearth.settings.isgrd==0,
					error('solidearthsettings checkconsistency error message; if bottom pressure grd patterns are requested, solidearth settings class should have isgrd flag on');
			end

		end % }}}
		function marshall(self,prefix,md,fid) % {{{
			WriteData(fid,prefix,'object',self,'fieldname','reltol','name','md.solidearth.settings.reltol','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','abstol','name','md.solidearth.settings.abstol','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','maxiter','name','md.solidearth.settings.maxiter','format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','selfattraction','name','md.solidearth.settings.selfattraction','format','Boolean');
			WriteData(fid,prefix,'object',self,'fieldname','elastic','name','md.solidearth.settings.elastic','format','Boolean');
			WriteData(fid,prefix,'object',self,'fieldname','viscous','name','md.solidearth.settings.viscous','format','Boolean');
			WriteData(fid,prefix,'object',self,'fieldname','rotation','name','md.solidearth.settings.rotation','format','Boolean');
			WriteData(fid,prefix,'object',self,'fieldname','grdocean','name','md.solidearth.settings.grdocean','format','Boolean');
			WriteData(fid,prefix,'object',self,'fieldname','ocean_area_scaling','name','md.solidearth.settings.ocean_area_scaling','format','Boolean');
			WriteData(fid,prefix,'object',self,'fieldname','runfrequency','name','md.solidearth.settings.runfrequency','format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','degacc','name','md.solidearth.settings.degacc','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','timeacc','name','md.solidearth.settings.timeacc','format','Double','scale',md.constants.yts);
			WriteData(fid,prefix,'object',self,'fieldname','horiz','name','md.solidearth.settings.horiz','format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','sealevelloading','name','md.solidearth.settings.sealevelloading','format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','isgrd','name','md.solidearth.settings.isgrd','format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','compute_bp_grd','name','md.solidearth.settings.compute_bp_grd','format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','grdmodel','name','md.solidearth.settings.grdmodel','format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','cross_section_shape','name','md.solidearth.settings.cross_section_shape','format','Integer');
		end % }}}
		function self = extrude(self,md) % {{{
		end % }}}
		function savemodeljs(self,fid,modelname) % {{{
			% TODO: Update the following fields so that coverage is complete
			writejsdouble(fid,[modelname '.solidearth.settings.maxiter'],self.maxiter);
			writejsdouble(fid,[modelname '.solidearth.settings.reltol'],self.reltol);
			writejsdouble(fid,[modelname '.solidearth.settings.abstol'],self.abstol);
			writejsdouble(fid,[modelname '.solidearth.settings.selfattraction'],self.selfattraction);
			writejsdouble(fid,[modelname '.solidearth.settings.elastic'],self.elastic);
			writejsdouble(fid,[modelname '.solidearth.settings.viscous'],self.viscous);
			writejsdouble(fid,[modelname '.solidearth.settings.rotation'],self.rotation);
			writejsdouble(fid,[modelname '.solidearth.settings.grdocean'],self.grdocean);
			writejsdouble(fid,[modelname '.solidearth.settings.ocean_area_scaling'],self.ocean_area_scaling);
			writejsdouble(fid,[modelname '.solidearth.settings.run_frequency'],self.run_frequency);
			writejsdouble(fid,[modelname '.solidearth.settings.degacc'],self.degacc);
			writejsdouble(fid,[modelname '.solidearth.settings.timeacc'],self.timeacc);
			writejsdouble(fid,[modelname '.solidearth.settings.cross_section_shape'],self.cross_section_shape);
		end % }}}
	end
end
