%CALVINCREVASSEDEPTH class definition
%
%   Usage:
%      calvingcrevassedepth=calvingcrevassedepth();

classdef calvingcrevassedepth
	properties (SetAccess=public) 
		crevasse_opening_stress=1;
		crevasse_threshold     =1.;
		surface_hydrology_type = 1;
		water_height = 0.;
		hurst = 0.;
		sigma = 0;
		melt_supply =0;
	end
	methods
		function self = calvingcrevassedepth(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				case 1
					inputstruct=varargin{1};
					list1 = properties('calvingcrevassedepth');
					list2 = fieldnames(inputstruct);
					for i=1:length(list1)
						fieldname = list1{i};
						if ismember(fieldname,list2)
							self.(fieldname) = inputstruct.(fieldname);
						end
					end
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = extrude(self,md) % {{{
		end % }}}
		function self = setdefaultparameters(self) % {{{
			
			self.crevasse_threshold      = 1.;
			self.crevasse_opening_stress = 1;
			self.surface_hydrology_type = 0;
         self.water_height       = 0.;
			self.hurst = 0.41;
			self.sigma = 14.3;
			self.melt_supply = 0;
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{
			%Early return
			if (~strcmp(solution,'TransientSolution') | md.transient.ismovingfront==0), return; end

			md = checkfield(md,'fieldname','calving.crevasse_opening_stress','numel',[1],'values',[0,1,2]);
         		md = checkfield(md,'fieldname','calving.crevasse_threshold','numel',[1],'>',0,'<=',1);
			md = checkfield(md,'fieldname','calving.surface_hydrology_type','numel',[1],'values',[0,1,2]);
			if self.surface_hydrology_type == 0
				md = checkfield(md,'fieldname','calving.water_height','NaN',1,'Inf',1,'timeseries',1,'>=',0);
			elseif self.surface_hydrology_type ==1
				md = checkfield(md,'fieldname','calving.hurst','NaN',1,'size',[md.mesh.numberofvertices 1],'>=',0,'<=',1);
                                md = checkfield(md,'fieldname','calving.sigma','NaN',1,'size',[md.mesh.numberofvertices 1],'>=',0);
                                md = checkfield(md,'fieldname','calving.melt_supply','NaN',1,'Inf',1,'timeseries',1,'>=',0);
			elseif self.surface_hydrology_type ==2
				md = checkfield(md,'fieldname','calving.hurst','NaN',1,'size',[md.mesh.numberofvertices 1],'>=',0,'<=',1);
                                md = checkfield(md,'fieldname','calving.sigma','NaN',1,'size',[md.mesh.numberofvertices 1],'>=',0);
                                md = checkfield(md,'fieldname','calving.melt_supply','NaN',1,'Inf',1,'timeseries',1,'>=',0);
			end
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   Calving Pi parameters:'));
			fielddisplay(self,'crevasse_opening_stress','0: stress only in the ice-flow direction, 1: max principal, 2: buttressing based');
			fielddisplay(self,'crevasse_threshold','ratio of full thickness to calve (e.g. 0.75 is for 75% of the total ice thickness)');
			fielddisplay(self,'surface_hydrology_type','0: default crevasse depth calving, 1: use mean water depth param to prescibe water_height, 2: use mean lake depth param to prescribe as water_height')
			fielddisplay(self,'water_height','water height in the crevasse [m]');
			fielddisplay(self,'water_height','water height in the crevasse [m]');
                        fielddisplay(self,'hurst','Hurst exponent');
            		fielddisplay(self,'sigma','Roughness Amplitude [m]')
            		fielddisplay(self,'melt_supply','Amount of supply meltwater able to be distributed across the surface [m]')

		end % }}}
		function marshall(self,prefix,md,fid) % {{{
			yts=md.constants.yts;
			WriteData(fid,prefix,'name','md.calving.law','data',6,'format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','crevasse_opening_stress','format','Integer');
         		WriteData(fid,prefix,'object',self,'fieldname','crevasse_threshold','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','surface_hydrology_type','format','Integer');
			if self.surface_hydrology_type ==0
				WriteData(fid,prefix,'object',self,'fieldname','water_height','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			elseif self.surface_hydrology_type ==1
                            WriteData(fid,prefix,'object',self,'fieldname','hurst','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
                            WriteData(fid,prefix,'object',self,'fieldname','sigma','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
                            WriteData(fid,prefix,'object',self,'fieldname','melt_supply','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			elseif self.surface_hydrology_type ==2
                                WriteData(fid,prefix,'object',self,'fieldname','hurst','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
                            WriteData(fid,prefix,'object',self,'fieldname','sigma','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
                            WriteData(fid,prefix,'object',self,'fieldname','melt_supply','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
                        end
		end % }}}
	end
end
