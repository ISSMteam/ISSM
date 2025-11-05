%SMBpddGCM Class definition
%
%   Usage:
%      SMBpddGCM=SMBpddGCM();

classdef SMBpddGCM
	properties (SetAccess=public)

		x_grid                = NaN;
		y_grid                = NaN;
		precipitation         = NaN;
		temperature           = NaN;
		enhance_factor        = NaN;
		lapserates            = NaN;
		allsolidtemperature   = 0+273.15;
		allliquidtemperature  = 2+273.15;
		ddf_snow              = 4.5/1000/24/3600; 
		ddf_ice               = 6.5/1000/24/3600;
      steps_per_step        = 1;
      averaging             = 0;
		ref_surf					 = NaN;
		requested_outputs     = {};
	end
	methods
		function self = SMBpddGCM(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = extrude(self,md) % {{{
			self.enhance_factor=project3d(md,'vector',self.enhance_factor,'type','node');
			self.lapserates=project3d(md,'vector',self.lapserates,'type','node');

		end % }}}
		function list = defaultoutputs(self,md) % {{{
			list = {'SmbMassBalance'};
		end % }}}
		function self = initialize(self,md) % {{{

			if isnan(self.enhance_factor),
				self.enhance_factor=zeros(md.mesh.numberofvertices,1);
				disp('      no SMBpddGCM.enhance_factor specified: values set as zero');
			end
			if isnan(self.lapserates),
				self.lapserates=0.0065*ones(md.mesh.numberofvertices,1);
				disp('      no SMBpddGCM.lapserates specified: values set as 0.0065 °C/m');
			end
			if isnan(self.ref_surf),
				self.ref_surf = 0*ones(md.mesh.numberofvertices,1);
				disp('      no SMBpddGCM.ref_surf specified: values set as sea-level 0 m');
			end
		end % }}}
		function self = setdefaultparameters(self) % {{{

			self.requested_outputs={'default'};

		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			if (strcmp(solution,'TransientSolution') & md.transient.issmb == 0), return; end

			if ismember('MasstransportAnalysis',analyses),
				md = checkfield(md,'fieldname','smb.x_grid','NaN',1,'Inf',1,'size',[NaN 1]);
				md = checkfield(md,'fieldname','smb.y_grid','NaN',1,'Inf',1,'size',[NaN 1]);
				Nx = numel(self.x_grid);
				Ny = numel(self.y_grid);
				md = checkfield(md,'fieldname','smb.precipitation','NaN',1,'Inf',1,'size',[Nx*Ny+1, NaN]);
				md = checkfield(md,'fieldname','smb.temperature','NaN',1,'Inf',1,'size',[Nx*Ny+1, NaN]);
				md = checkfield(md,'fieldname','smb.lapserates','>=',0,'NaN',1,'Inf',1,'size',[md.mesh.numberofvertices 1]);
				md = checkfield(md,'fieldname','smb.enhance_factor','>=',0,'NaN',1,'Inf',1,'size',[md.mesh.numberofvertices 1]);
				md = checkfield(md,'fieldname','smb.ref_surf','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices 1]);
				md = checkfield(md,'fieldname','smb.allsolidtemperature','NaN',1,'Inf',1);
				md = checkfield(md,'fieldname','smb.allliquidtemperature','NaN',1,'Inf',1);
				md = checkfield(md,'fieldname','smb.ddf_snow','>=',0,'NaN',1,'Inf',1);
				md = checkfield(md,'fieldname','smb.ddf_ice','>=',0,'NaN',1,'Inf',1);
			end
			md = checkfield(md,'fieldname','smb.steps_per_step','>=',1,'numel',[1]);
         md = checkfield(md,'fieldname','smb.averaging','numel',[1],'values',[0 1 2]);
			md = checkfield(md,'fieldname','smb.requested_outputs','stringrow',1);

		end % }}}
		function disp(self) % {{{
			disp(sprintf('   surface forcings parameters:'));

			disp(sprintf('\n   PDD Model based on downscaling of GCM climate data :'));
			fielddisplay(self,'x_grid','x coordinates of the GCM climate grid data');
			fielddisplay(self,'y_grid','y coordinates of the GCM climate grid data');
			fielddisplay(self,'precipitation', 'monthly surface precipitation ((nx*ny+1)*t) [m/yr water eq]');
			fielddisplay(self,'temperature', 'surface temperature ((nx*ny+1)*t) [K]');
			fielddisplay(self,'enhance_factor', 'melting enhance factor of degree-day factor for debris');
			fielddisplay(self,'lapserates', 'lapse rate [K/m]');
			fielddisplay(self,'allsolidtemperature', 'temperature where all precipitation is solid [K]');
			fielddisplay(self,'allliquidtemperature', 'temperature where all precipitation is liquid [K]');
			fielddisplay(self,'ddf_snow', 'DDF for snow (m w.e./s/K), Litrature: 4.1 ± 1.5 mm/d/K');
			fielddisplay(self,'ddf_ice', 'DDF for ice (m w.e./s/K), Litrature: 8.0 ± 3.4 mm/d/K');
         fielddisplay(self, 'steps_per_step', 'number of smb steps per time step');
         fielddisplay(self, 'averaging', 'averaging methods from short to long steps');
         fielddisplay(self, 'ref_surf', 'reference surface elevation for downsampling the GCM temperature');
         disp(sprintf('%51s  0: Arithmetic (default)',' '));
         disp(sprintf('%51s  1: Geometric',' '));
         disp(sprintf('%51s  2: Harmonic',' '));
			fielddisplay(self,'requested_outputs','additional outputs requested (TemperaturePDD, SmbAccumulation, SmbMelt)');
		end % }}}
		function marshall(self,prefix,md,fid) % {{{

			yts=md.constants.yts;
			Nx = numel(self.x_grid);
			Ny = numel(self.y_grid);

			% everything in the bin is in SI
			if isnan(self.ref_surf),
				% prepare ref surface for the GCM grid, using the initial md.geometry.surface
				grid_surf = InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,md.geometry.surface,self.x_grid,self.y_grid,NaN);
				self.ref_surf = InterpFromGridToMesh(self.x_grid,self.y_grid,grid_surf,md.mesh.x,md.mesh.y,0);
			end

			WriteData(fid,prefix,'name','md.smb.model','data',15,'format','Integer');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','x_grid','format','DoubleMat','mattype',3);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','y_grid','format','DoubleMat','mattype',3);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','precipitation','format','DoubleMat','mattype',3,'timeserieslength', Nx*Ny+1, 'scale',1./yts, 'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','temperature','format','DoubleMat','mattype',3,'timeserieslength', Nx*Ny+1, 'scale',1., 'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','enhance_factor','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','lapserates','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','allsolidtemperature','format','Double');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','allliquidtemperature','format','Double');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','ddf_snow','format','Double');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','ddf_ice','format','Double');

			WriteData(fid, prefix, 'object', self, 'fieldname', 'steps_per_step', 'format', 'Integer');
         WriteData(fid, prefix, 'object', self, 'fieldname', 'averaging', 'format', 'Integer')
			%process requested outputs
			outputs = self.requested_outputs;
			pos  = find(ismember(outputs,'default'));
			if ~isempty(pos),
				outputs(pos) = [];                         %remove 'default' from outputs
				outputs      = [outputs defaultoutputs(self,md)]; %add defaults
			end
			WriteData(fid,prefix,'data',outputs,'name','md.smb.requested_outputs','format','StringArray');

			WriteData(fid,prefix,'object',self,'class','smb','fieldname','ref_surf','format','DoubleMat','mattype',1);
		end % }}}
	end
end
