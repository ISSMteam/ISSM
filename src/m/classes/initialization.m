%INITIALIZATION class definition
%
%   Usage:
%      initialization=initialization();

classdef initialization
	properties (SetAccess=public) 
		vx                  = NaN;
		vy                  = NaN;
		vz                  = NaN;
		vel                 = NaN;
		pressure            = NaN;
		temperature         = NaN;
		enthalpy            = NaN;
		waterfraction       = NaN;
		sediment_head       = NaN;
		epl_head            = NaN;
		epl_thickness       = NaN;
		watercolumn         = NaN;
		hydraulic_potential = NaN;
		channelarea         = NaN;
		sealevel            = NaN;
		bottompressure      = NaN;
		dsl                 = NaN;
		str                 = NaN;
		sample              = NaN;
		debris              = NaN;
		age                 = NaN;
	end
	methods
		function self = initialization(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = setdefaultparameters(self) % {{{
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{
			if ismember('StressbalanceAnalysis',analyses) & ~(strcmp(solution,'TransientSolution') & md.transient.isstressbalance == 0),
				if numel(md.initialization.vx)>1 | numel(md.initialization.vy)>1
					md = checkfield(md,'fieldname','initialization.vx','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices 1]);
					md = checkfield(md,'fieldname','initialization.vy','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices 1]);
				end
			end
			if ismember('MasstransportAnalysis',analyses) & ~(strcmp(solution,'TransientSolution') & md.transient.ismasstransport == 0),
				md = checkfield(md,'fieldname','initialization.vx','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices 1]);
				md = checkfield(md,'fieldname','initialization.vy','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices 1]);
			end
			if ismember('OceantransportAnalysis',analyses) 
				if strcmp(solution,'TransientSolution') & md.transient.isslc & md.transient.isoceantransport,
					md = checkfield(md,'fieldname','initialization.bottompressure','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices 1]);
					md = checkfield(md,'fieldname','initialization.dsl','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices 1]);
					md = checkfield(md,'fieldname','initialization.str','NaN',1,'Inf',1,'size',[1 1]);
				end
			end
			if ismember('BalancethicknessAnalysis',analyses) & strcmp(solution,'BalancethicknessSolution'),
				md = checkfield(md,'fieldname','initialization.vx','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices 1]);
				md = checkfield(md,'fieldname','initialization.vy','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices 1]);
				%Triangle with zero velocity
				if any(sum(abs(md.initialization.vx(md.mesh.elements)),2)==0 & sum(abs(md.initialization.vy(md.mesh.elements)),2)==0 & min(md.mask.ice_levelset(md.mesh.elements),[],2)<0)
					md = checkmessage(md,'at least one triangle has all its vertices with a zero velocity');
				end
			end
			if ismember('ThermalAnalysis',analyses) & ~(strcmp(solution,'TransientSolution') & md.transient.isthermal == 0),
				md = checkfield(md,'fieldname','initialization.vx','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices 1]);
				md = checkfield(md,'fieldname','initialization.vy','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices 1]);
				if dimension(md.mesh)==3
					md = checkfield(md,'fieldname','initialization.vz','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices 1]);
				end
				md = checkfield(md,'fieldname','initialization.pressure','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices 1]);
				md = checkfield(md,'fieldname','initialization.temperature','NaN',1,'Inf',1,'size','universal');
			end
			if (ismember('EnthalpyAnalysis',analyses) & md.thermal.isenthalpy)
				md = checkfield(md,'fieldname','initialization.waterfraction','>=',0,'size','universal');
				md = checkfield(md,'fieldname','initialization.watercolumn'  ,'>=',0,'size','universal');
				pos=find(md.initialization.waterfraction>0.);
				if(~isempty(pos)),
					%md = checkfield(md,'fieldname', 'delta Tpmp', 'field',abs(md.initialization.temperature(pos)-(md.materials.meltingpoint-md.materials.beta*md.initialization.pressure(pos))),'<',1e-11,...
					%'message','set temperature to pressure melting point at locations with waterfraction>0');
				end
			end
			if ismember('HydrologyShreveAnalysis',analyses),
				if isa(md.hydrology,'hydrologyshreve'),
					if (strcmp(solution,'TransientSolution') & md.transient.ishydrology) | strcmp(solution,'HydrologySolution'),
						md = checkfield(md,'fieldname','initialization.watercolumn','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices 1]);
					end
				end
			end
			if ismember('HydrologyTwsAnalysis',analyses),
				if isa(md.hydrology,'hydrologytws'),
					md = checkfield(md,'fieldname','initialization.watercolumn','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices 1]);
				end
			end
			if ismember('SealevelchangeAnalysis',analyses),
				if strcmp(solution,'TransientSolution') & md.transient.isslc,
					md = checkfield(md,'fieldname','initialization.sealevel','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices 1]);
				end
			end
			if ismember('HydrologyGlaDSAnalysis',analyses),
				if isa(md.hydrology,'hydrologyglads'),
					md = checkfield(md,'fieldname','initialization.watercolumn','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices 1]);
					md = checkfield(md,'fieldname','initialization.hydraulic_potential','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices 1]);
					md = checkfield(md,'fieldname','initialization.channelarea','NaN',1,'Inf',1,'>=',0,'size',[md.mesh.numberofedges 1]);
				end
			end
			if ismember('HydrologyDCInefficientAnalysis',analyses),
				if isa(md.hydrology,'hydrologydc'),
					md = checkfield(md,'fieldname','initialization.sediment_head','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices 1]);
				end
			end
			if ismember('HydrologyDCEfficientAnalysis',analyses),
				if isa(md.hydrology,'hydrologydc'),
					if md.hydrology.isefficientlayer==1,
						md = checkfield(md,'fieldname','initialization.epl_head','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices 1]);
						md = checkfield(md,'fieldname','initialization.epl_thickness','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices 1]);
					end
				end
			end
			if ismember('SamplingAnalysis',analyses) & ~(strcmp(solution,'TransientSolution') & md.transient.issampling == 0),
				if ~isnan(md.initialization.sample)
					md = checkfield(md,'fieldname','initialization.sample','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices 1]);
				end
			end
			if ismember('DebrisAnalysis',analyses),
				if ~isnan(md.initialization.debris)
					md = checkfield(md,'fieldname','initialization.debris','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices 1]);
				end
			end
			if ismember('AgeAnalysis',analyses),
				if ~isnan(md.initialization.age)
					md = checkfield(md,'fieldname','initialization.age','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices 1]);
				end
			end
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   initial field values:'));

			fielddisplay(self,'vx','x component of velocity [m/yr]');
			fielddisplay(self,'vy','y component of velocity [m/yr]');
			fielddisplay(self,'vz','z component of velocity [m/yr]');
			fielddisplay(self,'vel','velocity norm [m/yr]');
			fielddisplay(self,'pressure','pressure field [Pa]');
			fielddisplay(self,'temperature','temperature [K]');
			fielddisplay(self,'enthalpy','enthalpy [J]');
			fielddisplay(self,'waterfraction','fraction of water in the ice');
			fielddisplay(self,'sediment_head','sediment water head of subglacial system [m]');
			fielddisplay(self,'epl_head','epl water head of subglacial system [m]');
			fielddisplay(self,'epl_thickness','epl layer thickness [m]');
			fielddisplay(self,'watercolumn','subglacial water sheet thickness (for Shreve and GlaDS) [m]');
			fielddisplay(self,'hydraulic_potential','Hydraulic potential (for GlaDS) [Pa]');
			fielddisplay(self,'channelarea','subglacial water channel area (for GlaDS) [m2]');
			fielddisplay(self,'sample','Realization of a Gaussian random field');
			fielddisplay(self,'bottompressure','Bottom pressures');
			fielddisplay(self,'dsl','Dynamic sea level.');
			fielddisplay(self,'str','Steric sea level.');
			fielddisplay(self,'debris','Surface debris layer [m]');
			fielddisplay(self,'age','Initial age [yr]');
		end % }}}
		function marshall(self,prefix,md,fid) % {{{

			yts=md.constants.yts;

			WriteData(fid,prefix,'object',self,'fieldname','vx','format','DoubleMat','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1,'yts',yts);
			WriteData(fid,prefix,'object',self,'fieldname','vy','format','DoubleMat','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1,'yts',yts);
			WriteData(fid,prefix,'object',self,'fieldname','vz','format','DoubleMat','mattype',1,'scale',1./yts);
			WriteData(fid,prefix,'object',self,'fieldname','pressure','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'fieldname','sealevel','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',yts);
			WriteData(fid,prefix,'object',self,'fieldname','bottompressure','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'fieldname','str','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'fieldname','dsl','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'fieldname','temperature','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'fieldname','waterfraction','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'fieldname','sediment_head','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'fieldname','epl_head','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'fieldname','epl_thickness','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'fieldname','watercolumn','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'fieldname','channelarea','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'fieldname','hydraulic_potential','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'fieldname','sample','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'fieldname','debris','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'fieldname','age','format','DoubleMat','mattype',1,'scale',yts);

			if md.thermal.isenthalpy,
				if numel(self.enthalpy) <= 1,
					%reconstruct enthalpy
					tpmp = md.materials.meltingpoint - md.materials.beta*md.initialization.pressure;
					pos  = find(md.initialization.waterfraction>0.);
					enthalpy      = md.materials.heatcapacity*(md.initialization.temperature-md.constants.referencetemperature);
					enthalpy(pos) = md.materials.heatcapacity*(tpmp(pos) - md.constants.referencetemperature) + md.materials.latentheat*md.initialization.waterfraction(pos);
				else
					enthalpy = self.enthalpy;
				end
				WriteData(fid,prefix,'data',enthalpy,'format','DoubleMat','mattype',1,'name','md.initialization.enthalpy');
			end
		end % }}}
		function self = extrude(self,md) % {{{
			self.vx=project3d(md,'vector',self.vx,'type','node');
			self.vy=project3d(md,'vector',self.vy,'type','node');
			self.vz=project3d(md,'vector',self.vz,'type','node');
			self.vel=project3d(md,'vector',self.vel,'type','node');
			self.temperature=project3d(md,'vector',self.temperature,'type','node');
			self.enthalpy=project3d(md,'vector',self.enthalpy,'type','node');
			self.waterfraction=project3d(md,'vector',self.waterfraction,'type','node');
			self.watercolumn=project3d(md,'vector',self.watercolumn,'type','node','layer',1);
			self.sediment_head=project3d(md,'vector',self.sediment_head,'type','node','layer',1);
			self.epl_head=project3d(md,'vector',self.epl_head,'type','node','layer',1);
			self.epl_thickness=project3d(md,'vector',self.epl_thickness,'type','node','layer',1);
			self.sealevel=project3d(md,'vector',self.sealevel,'type','node','layer',1);
			self.bottompressure=project3d(md,'vector',self.bottompressure,'type','node','layer',1);
			self.dsl=project3d(md,'vector',self.dsl,'type','node','layer',1);
			self.str=project3d(md,'vector',self.str,'type','node','layer',1);
			self.str=project3d(md,'vector',self.debris,'type','node','layer',1);
			self.str=project3d(md,'vector',self.age,'type','node','layer',1);

			%Lithostatic pressure by default
			self.pressure=md.constants.g*md.materials.rho_ice*(md.geometry.surface-md.mesh.z);
		end % }}}
		function savemodeljs(self,fid,modelname) % {{{

			writejs1Darray(fid,[modelname '.initialization.vx'],self.vx);
			writejs1Darray(fid,[modelname '.initialization.vy'],self.vy);
			writejs1Darray(fid,[modelname '.initialization.vz'],self.vz);
			writejs1Darray(fid,[modelname '.initialization.vel'],self.vel);
			writejs1Darray(fid,[modelname '.initialization.pressure'],self.pressure);
			writejs1Darray(fid,[modelname '.initialization.temperature'],self.temperature);
			writejs1Darray(fid,[modelname '.initialization.enthalpy'],self.enthalpy);
			writejs1Darray(fid,[modelname '.initialization.waterfraction'],self.waterfraction);
			writejs1Darray(fid,[modelname '.initialization.sediment_head'],self.sediment_head);
			writejs1Darray(fid,[modelname '.initialization.epl_head'],self.epl_head);
			writejs1Darray(fid,[modelname '.initialization.epl_thickness'],self.epl_thickness);
			writejs1Darray(fid,[modelname '.initialization.watercolumn'],self.watercolumn);
			writejs1Darray(fid,[modelname '.initialization.hydraulic_potential'],self.hydraulic_potential);
			writejs1Darray(fid,[modelname '.initialization.channel'],self.channelarea);
			writejs1Darray(fid,[modelname '.initialization.sample'],self.sample);
			writejs1Darray(fid,[modelname '.initialization.debris'],self.debris);
			writejs1Darray(fid,[modelname '.initialization.age'],self.age);

		end % }}}
	end
end
