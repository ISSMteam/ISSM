%SMBd18opdd Class definition
%
%   Usage:
%      SMBd18opdd=SMBd18opdd();

classdef SMBd18opdd
	properties (SetAccess=public)

		desfac                    = 0;
		s0p                       = NaN;
		s0t                       = NaN;
		rlaps                     = 0;
		rlapslgm                  = 0;
		dpermil                   = 0;
		f                         = 0;
		Tdiff                     = NaN;
		sealev                    = NaN;
		ismungsm                  = 0;
		isd18opd                  = 0;
		issetpddfac               = 0;
		istemperaturescaled       = 1;
		isprecipscaled            = 1;
		delta18o                  = NaN;
		delta18o_surface          = NaN;
		temperatures_presentday   = NaN;
		precipitations_presentday = NaN;
		temperatures_reconstructed = NaN;
		precipitations_reconstructed = NaN;
		pddfac_snow               = NaN;
		pddfac_ice                = NaN;
		steps_per_step            = 1;
		averaging                 = 0;
		requested_outputs         = {};
	end
	methods
		function self = SMBd18opdd(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = extrude(self,md) % {{{
			if(self.isd18opd),self.temperatures_presentday=project3d(md,'vector',self.temperatures_presentday,'type','node');end
			if(self.isd18opd),self.precipitations_presentday=project3d(md,'vector',self.precipitations_presentday,'type','node');end
			if(self.istemperaturescaled==0),self.temperatures_reconstructed=project3d(md,'vector',self.temperatures_reconstructed,'type','node');end
			if(self.isprecipscaled==0),self.precipitations_reconstructed=project3d(md,'vector',self.precipitations_reconstructed,'type','node');end
			self.s0p=project3d(md,'vector',self.s0p,'type','node');
			self.s0t=project3d(md,'vector',self.s0t,'type','node');

		end % }}}
		function list = defaultoutputs(self,md) % {{{

			list = {'SmbMassBalance'};

		end % }}}
		function self = initialize(self,md) % {{{

			if isnan(self.s0p),
				self.s0p=zeros(md.mesh.numberofvertices,1);
				disp('      no SMBd18opdd.s0p specified: values set as zero');
			end
			if isnan(self.s0t),
				self.s0t=zeros(md.mesh.numberofvertices,1);
				disp('      no SMBd18opdd.s0t specified: values set as zero');
			end

		end % }}}
		function self = setdefaultparameters(self) % {{{

			self.ismungsm   = 0;
			self.isd18opd   = 1;
			self.istemperaturescaled = 1;
			self.isprecipscaled = 1;
			self.desfac     = 0.5;
			self.rlaps      = 6.5;
			self.rlapslgm   = 6.5;
			self.dpermil    = 2.4;
			self.f          = 0.169;
			self.issetpddfac = 0;
			self.requested_outputs={'default'};
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			if ismember('MasstransportAnalysis',analyses),
				md = checkfield(md,'fieldname','smb.desfac','<=',1,'numel',1);
				md = checkfield(md,'fieldname','smb.s0p','>=',0,'NaN',1,'Inf',1,'size',[md.mesh.numberofvertices 1]);
				md = checkfield(md,'fieldname','smb.s0t','>=',0,'NaN',1,'Inf',1,'size',[md.mesh.numberofvertices 1]);
				md = checkfield(md,'fieldname','smb.rlaps','>=',0,'numel',1);
				md = checkfield(md,'fieldname','smb.rlapslgm','>=',0,'numel',1);

				if(self.isd18opd==1)
					md = checkfield(md,'fieldname','smb.temperatures_presentday','size',[md.mesh.numberofvertices 12],'NaN',1,'Inf',1);
					md = checkfield(md,'fieldname','smb.precipitations_presentday','size',[md.mesh.numberofvertices 12],'NaN',1,'Inf',1);
					md = checkfield(md,'fieldname','smb.delta18o','NaN',1,'Inf',1,'size',[2,NaN],'singletimeseries',1);
					md = checkfield(md,'fieldname','smb.dpermil','>=',0,'numel',1);
					md = checkfield(md,'fieldname','smb.f','>=',0,'numel',1);
					if(self.istemperaturescaled==0)
						lent=size(self.temperatures_reconstructed,2);
						multt=ceil(lent/12)*12;
						md = checkfield(md,'fieldname','smb.temperatures_reconstructed','size',[md.mesh.numberofvertices+1 multt],'NaN',1,'Inf',1,'timeseries',1);
					end
					if(self.isprecipscaled==0)
						lenp=size(self.temperatures_reconstructed,2);
						multp=ceil(lent/12)*12;
						md = checkfield(md,'fieldname','smb.precipitations_reconstructed','size',[md.mesh.numberofvertices+1 multt],'NaN',1,'Inf',1,'timeseries',1);
					end
				end
				if(self.issetpddfac==1)
					md = checkfield(md,'fieldname','smb.pddfac_snow','>=',0,'NaN',1,'Inf',1);
					md = checkfield(md,'fieldname','smb.pddfac_ice','>=',0,'NaN',1,'Inf',1);
				end
			end
			md = checkfield(md,'fieldname','smb.steps_per_step','>=',1,'numel',[1]);
			md = checkfield(md,'fieldname','smb.averaging','numel',[1],'values',[0 1 2]);
			md = checkfield(md,'fieldname','smb.requested_outputs','stringrow',1);
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   surface forcings parameters:'));

			disp(sprintf('\n   PDD and deltaO18 parameters:'));
			fielddisplay(self,'isd18opd','is delta18o parametrisation from present day temperature and precipitation activated (0 or 1, default is 0)');
			fielddisplay(self,'issetpddfac','is user passing in defined pdd factors (0 or 1, default is 0)');
			fielddisplay(self,'desfac','desertification elevation factor (between 0 and 1, default is 0.5) [m]');
			fielddisplay(self,'s0p','should be set to elevation from precip source (between 0 and a few 1000s m, default is 0) [m]');
			fielddisplay(self,'s0t','should be set to elevation from temperature source (between 0 and a few 1000s m, default is 0) [m]');
			fielddisplay(self,'rlaps','present day lapse rate [degree/km]');
			if(self.isd18opd==1)
				fielddisplay(self,'temperatures_presentday','monthly present day surface temperatures [K], required if delta18o/mungsm/d18opd is activated');
				fielddisplay(self,'precipitations_presentday','monthly surface precipitation [m/yr water eq], required if delta18o/mungsm/d18opd is activated');
				fielddisplay(self,'istemperaturescaled','if delta18o parametrisation from present day temperature and precipitation is activated, is temperature scaled to delta18o value? (0 or 1, default is 1)');
				fielddisplay(self,'isprecipscaled','if delta18o parametrisation from present day temperature and precipitation is activated, is precip scaled to delta18o value? (0 or 1, default is 1)');
				if(self.istemperaturescaled==0)
					fielddisplay(self,'temperatures_reconstructed','monthly historical surface temperatures [K], required if delta18o/mungsm/d18opd is activated and istemperaturescaled is not activated');
				end
				if(self.isprecipscaled==0)
					fielddisplay(self,'precipitations_reconstructed','monthly historical precipitation [m/yr water eq], required if delta18o/mungsm/d18opd is activated and isprecipscaled is not activated');
				end
				fielddisplay(self,'delta18o','delta18o [per mil], required if pdd is activated and d18opd activated');
				fielddisplay(self,'dpermil','degree per mil, required if d18opd is activated');
				fielddisplay(self,'f','precip/temperature scaling factor, required if d18opd is activated');
			end
			if(self.issetpddfac==1)
				fielddisplay(self,'pddfac_snow','Pdd factor for snow for all the domain [mm ice equiv/day/degree C]');
				fielddisplay(self,'pddfac_ice','Pdd factor for ice for all the domain [mm ice equiv/day/degree C]');
			end
			fielddisplay(self, 'steps_per_step', 'number of smb steps per time step');
			fielddisplay(self, 'averaging', 'averaging methods from short to long steps');
			disp(sprintf('%51s  0: Arithmetic (default)',' '));
			disp(sprintf('%51s  1: Geometric',' '));
			disp(sprintf('%51s  2: Harmonic',' '));
			fielddisplay(self,'requested_outputs','additional outputs requested');
			% No need to display rlapslgm, Tdiff, ismungsm
			% as they are not used in this case but are still needed as default values in
			% PositiveDegreeDay (Tria.cpp) used in that case
		end % }}}
		function marshall(self,prefix,md,fid) % {{{

			yts=md.constants.yts;

			WriteData(fid,prefix,'name','md.smb.model','data',5,'format','Integer');

			WriteData(fid,prefix,'object',self,'class','smb','fieldname','ismungsm','format','Boolean');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','isd18opd','format','Boolean');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','issetpddfac','format','Boolean');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','desfac','format','Double');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','s0p','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','s0t','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','rlaps','format','Double');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','rlapslgm','format','Double');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','Tdiff','format','DoubleMat','mattype',1,'timeserieslength',2,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','sealev','format','DoubleMat','mattype',1,'timeserieslength',2,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'fieldname','steps_per_step','format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','averaging','format','Integer');

			if self.isd18opd
				WriteData(fid,prefix,'object',self,'class','smb','fieldname','temperatures_presentday','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
				WriteData(fid,prefix,'object',self,'class','smb','fieldname','precipitations_presentday','format','DoubleMat','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
				WriteData(fid,prefix,'object',self,'class','smb','fieldname','istemperaturescaled','format','Boolean');
				WriteData(fid,prefix,'object',self,'class','smb','fieldname','isprecipscaled','format','Boolean');
				WriteData(fid,prefix,'object',self,'class','smb','fieldname','delta18o','format','DoubleMat','mattype',1,'timeserieslength',2,'yts',md.constants.yts);
				WriteData(fid,prefix,'object',self,'class','smb','fieldname','dpermil','format','Double');
				WriteData(fid,prefix,'object',self,'class','smb','fieldname','f','format','Double');
				if self.istemperaturescaled==0
					WriteData(fid,prefix,'object',self,'class','smb','fieldname','temperatures_reconstructed','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
				end
				if self.isprecipscaled==0
					WriteData(fid,prefix,'object',self,'class','smb','fieldname','precipitations_reconstructed','format','DoubleMat','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
				end
			end
			if self.issetpddfac==1
				WriteData(fid,prefix,'object',self,'class','smb','fieldname','pddfac_snow','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
				WriteData(fid,prefix,'object',self,'class','smb','fieldname','pddfac_ice','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			end

			%process requested outputs
			outputs = self.requested_outputs;
			pos  = find(ismember(outputs,'default'));
			if ~isempty(pos),
				outputs(pos) = [];                         %remove 'default' from outputs
				outputs      = [outputs defaultoutputs(self,md)]; %add defaults
			end
			WriteData(fid,prefix,'data',outputs,'name','md.smb.requested_outputs','format','StringArray');

		end % }}}
	end
end
