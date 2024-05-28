%PLUME BASAL FORCINGS class definition
%
%   Usage:
%      plumebasalforcings=plumebasalforcings();

classdef plumebasalforcings
	properties (SetAccess=public) 
		floatingice_melting_rate  = NaN;
		groundedice_melting_rate  = NaN;
		mantleconductivity        = NaN;
		nusselt                   = NaN;
		dtbg                      = NaN;
		plumeradius               = NaN;
		topplumedepth             = NaN;
		bottomplumedepth          = NaN;
		plumex                    = NaN;
		plumey                    = NaN;
		crustthickness            = NaN;
		uppercrustthickness       = NaN;
		uppercrustheat            = NaN;
		lowercrustheat            = NaN;

	end
	methods
		function self = extrude(self,md) % {{{
			self.groundedice_melting_rate=project3d(md,'vector',self.groundedice_melting_rate,'type','node','layer',1); 
			self.floatingice_melting_rate=project3d(md,'vector',self.floatingice_melting_rate,'type','node','layer',1); 
		end % }}}
		function self = plumebasalforcings(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				case 1
					self=structtoobj(plumebasalforcings(),varargin{1});
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = initialize(self,md) % {{{

			if isnan(self.groundedice_melting_rate),
				self.groundedice_melting_rate=zeros(md.mesh.numberofvertices,1);
				disp('      no basalforcings.groundedice_melting_rate specified: values set as zero');
			end
			if isnan(self.floatingice_melting_rate),
				self.floatingice_melting_rate=zeros(md.mesh.numberofvertices,1);
				disp('      no basalforcings.floatingice_melting_rate specified: values set as zero');
			end

		end % }}}
		function self = setdefaultparameters(self) % {{{

			%default values for melting parameterization
			self.mantleconductivity     = 2.2;
			self.nusselt                = 300;
			self.dtbg                   = 11/1000.;
			self.plumeradius            = 100000;
			self.topplumedepth          = 10000;
			self.bottomplumedepth       = 1050000;
			self.crustthickness         = 30000;
			self.uppercrustthickness    = 14000;
			self.uppercrustheat         = 1.7*10^-6;
			self.lowercrustheat         = 0.4*10^-6;

		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			if ismember('MasstransportAnalysis',analyses) & ~(strcmp(solution,'TransientSolution') & md.transient.ismasstransport==0),
				md = checkfield(md,'fieldname','basalforcings.groundedice_melting_rate','NaN',1,'timeseries',1);
				md = checkfield(md,'fieldname','basalforcings.floatingice_melting_rate','NaN',1,'timeseries',1);
			end
			if ismember('BalancethicknessAnalysis',analyses),
				md = checkfield(md,'fieldname','basalforcings.groundedice_melting_rate','NaN',1,'size',[md.mesh.numberofvertices 1]);
				md = checkfield(md,'fieldname','basalforcings.floatingice_melting_rate','NaN',1,'size',[md.mesh.numberofvertices 1]);
			end
			if ismember('ThermalAnalysis',analyses) & ~(strcmp(solution,'TransientSolution') & md.transient.isthermal==0),
				md = checkfield(md,'fieldname','basalforcings.groundedice_melting_rate','NaN',1,'timeseries',1);
				md = checkfield(md,'fieldname','basalforcings.floatingice_melting_rate','NaN',1,'timeseries',1);
				md = checkfield(md,'fieldname','basalforcings.mantleconductivity','>=',0,'numel',1);
				md = checkfield(md,'fieldname','basalforcings.nusselt','>',0,'numel',1);
				md = checkfield(md,'fieldname','basalforcings.dtbg','>',0,'numel',1);
				md = checkfield(md,'fieldname','basalforcings.topplumedepth','>',0,'numel',1);
				md = checkfield(md,'fieldname','basalforcings.bottomplumedepth','>',0,'numel',1);
				md = checkfield(md,'fieldname','basalforcings.plumex','numel',1);
				md = checkfield(md,'fieldname','basalforcings.plumey','numel',1);
				md = checkfield(md,'fieldname','basalforcings.crustthickness','>',0,'numel',1);
				md = checkfield(md,'fieldname','basalforcings.uppercrustthickness','>',0,'numel',1);
				md = checkfield(md,'fieldname','basalforcings.uppercrustheat','>',0,'numel',1);
				md = checkfield(md,'fieldname','basalforcings.lowercrustheat','>',0,'numel',1);
			end
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   mantle plume basal melt parameterization:'));

			fielddisplay(self,'groundedice_melting_rate','basal melting rate (positive if melting) [m/yr]');
			fielddisplay(self,'floatingice_melting_rate','basal melting rate (positive if melting) [m/yr]');
			fielddisplay(self,'mantleconductivity','mantle heat conductivity [W/m^3]');
			fielddisplay(self,'nusselt','nusselt number, ratio of mantle to plume [1]');
			fielddisplay(self,'dtbg','background temperature gradient [degree/m]');
			fielddisplay(self,'plumeradius','radius of the mantle plume [m]');
			fielddisplay(self,'topplumedepth','depth of the mantle plume top below the crust [m]');
			fielddisplay(self,'bottomplumedepth','depth of the mantle plume base below the crust [m]');
			fielddisplay(self,'plumex','x coordinate of the center of the plume [m]');
			fielddisplay(self,'plumey','y coordinate of the center of the plume [m]');
			fielddisplay(self,'crustthickness','thickness of the crust [m]');
			fielddisplay(self,'uppercrustthickness','thickness of the upper crust [m]');
			fielddisplay(self,'uppercrustheat','volumic heat of the upper crust [w/m^3]');
			fielddisplay(self,'lowercrustheat','volumic heat of the lowercrust [w/m^3]');

		end % }}}
		function marshall(self,prefix,md,fid) % {{{

			yts=md.constants.yts;

			WriteData(fid,prefix,'name','md.basalforcings.model','data',4,'format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','floatingice_melting_rate','format','DoubleMat','name','md.basalforcings.floatingice_melting_rate','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts)
			WriteData(fid,prefix,'object',self,'fieldname','groundedice_melting_rate','format','DoubleMat','name','md.basalforcings.groundedice_melting_rate','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts)
			WriteData(fid,prefix,'object',self,'fieldname','mantleconductivity','format','Double')
			WriteData(fid,prefix,'object',self,'fieldname','nusselt','format','Double')
			WriteData(fid,prefix,'object',self,'fieldname','dtbg','format','Double')
			WriteData(fid,prefix,'object',self,'fieldname','plumeradius','format','Double')
			WriteData(fid,prefix,'object',self,'fieldname','topplumedepth','format','Double')
			WriteData(fid,prefix,'object',self,'fieldname','bottomplumedepth','format','Double')
			WriteData(fid,prefix,'object',self,'fieldname','plumex','format','Double')
			WriteData(fid,prefix,'object',self,'fieldname','plumey','format','Double')
			WriteData(fid,prefix,'object',self,'fieldname','crustthickness','format','Double')
			WriteData(fid,prefix,'object',self,'fieldname','uppercrustthickness','format','Double')
			WriteData(fid,prefix,'object',self,'fieldname','uppercrustheat','format','Double')
			WriteData(fid,prefix,'object',self,'fieldname','lowercrustheat','format','Double')
		end % }}}
	end
end
