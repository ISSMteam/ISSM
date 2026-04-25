%GROUNDINGLINE class definition
%
%   Usage:
%      groundingline=groundingline();

classdef groundingline
	properties (SetAccess=public) 
		migration              = '';
		friction_interpolation = '';
		melt_interpolation     = '';
		intrusion_distance     = 0;
		requested_outputs      = {};
	end
	methods
		function self = groundingline(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = setdefaultparameters(self) % {{{

			%Type of migration
			self.migration             = 'SubelementMigration';
			self.friction_interpolation= 'SubelementFriction1';
			self.melt_interpolation    = 'NoMeltOnPartiallyFloating';
			self.intrusion_distance    = 0;
			%default output
			self.requested_outputs     = {'default'};

		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			md = checkfield(md,'fieldname','groundingline.migration','values',{'None' 'SubelementMigration' 'AggressiveMigration' 'SoftMigration' 'Contact' 'GroundingOnly'});
			md = checkfield(md,'fieldname','groundingline.friction_interpolation','values',{'NoFrictionOnPartiallyFloating' 'SubelementFriction1' 'SubelementFriction2'});
			md = checkfield(md,'fieldname','groundingline.melt_interpolation','values',{'NoMeltOnPartiallyFloating' 'FullMeltOnPartiallyFloating' 'SubelementMelt1' 'SubelementMelt2' 'IntrusionMelt'});
			md = checkfield(md,'fieldname','groundingline.intrusion_distance','NaN',1,'Inf',1,'>=',0);
			md = checkfield(md,'fieldname','groundingline.requested_outputs','stringrow',1);

			if ~strcmp(self.migration,'None') & strcmp(solution,'TransientSolution') & md.transient.isgroundingline==1,
				if isnan(md.geometry.bed),
					md = checkmessage(md,['requesting grounding line migration, but bathymetry is absent!']);
				end
				pos=find(md.mask.ocean_levelset>0. & md.mask.ice_levelset<=0);
				if any(abs(md.geometry.base(pos)-md.geometry.bed(pos))>10^-10),
					md = checkmessage(md,['base not equal to bed on grounded ice!']);
				end
				pos=find(md.mask.ocean_levelset<=0. & md.mask.ice_levelset<=0);
				if any(md.geometry.bed(pos) - md.geometry.base(pos) > 10^-9),
					md = checkmessage(md,['bed superior to base on floating ice!']);
				end
			end

		end % }}}
		function list = defaultoutputs(self,md) % {{{
      
			list = {'Surface','Base','MaskOceanLevelset'};
         
      end % }}}
		function disp(self) % {{{
			disp(sprintf('   grounding line migration parameters:'));
			fielddisplay(self,'migration','type of grounding line migration: ''SoftMigration'',''SubelementMigration'',''AggressiveMigration'',''Contact'' or ''None''');
			fielddisplay(self,'friction_interpolation','type of friction interpolation for partially floating elements: ''NoFrictionOnPartiallyFloating'',''SubelementFriction1'', or ''SubelementFriction2''');
			fielddisplay(self,'melt_interpolation','type of melt interpolation for partially floating elements: ''NoMeltOnPartiallyFloating'',''FullMeltOnPartiallyFloating'',''SubelementMelt1'',''SubelementMelt2'' or ''IntrusionMelt''');
			fielddisplay(self,'intrusion_distance','distance of seawater intrusion from grounding line [m]');
			fielddisplay(self,'requested_outputs','additional outputs requested');

		end % }}}
		function marshall(self,prefix,md,fid) % {{{
			WriteData(fid,prefix,'data',self.migration,'name','md.groundingline.migration','format','String');
			WriteData(fid,prefix,'data',self.friction_interpolation,'name','md.groundingline.friction_interpolation','format','String');
			WriteData(fid,prefix,'data',self.melt_interpolation,'name','md.groundingline.melt_interpolation','format','String');
			WriteData(fid,prefix,'class','groundingline','object',self,'fieldname','intrusion_distance','format','DoubleMat','mattype',1);
			
			%process requested outputs
         outputs = self.requested_outputs;
         pos  = find(ismember(outputs,'default'));
         if ~isempty(pos),
            outputs(pos) = [];                         %remove 'default' from outputs
            outputs      = [outputs defaultoutputs(self,md)]; %add defaults
         end
			WriteData(fid,prefix,'data',outputs,'name','md.groundingline.requested_outputs','format','StringArray')
		end % }}}
		function savemodeljs(self,fid,modelname) % {{{
		
			writejsstring(fid,[modelname '.groundingline.migration'],self.migration);
			writejsstring(fid,[modelname '.groundingline.friction_interpolation'],self.friction_interpolation);
			writejsstring(fid,[modelname '.groundingline.melt_interpolation'],self.melt_interpolation);
			writejs1Darray(fid,[modelname '.groundingline.intrusion_distance'],self.intrusion_distance);
			writejscellstring(fid,[modelname '.groundingline.requested_outputs'],self.requested_outputs);

		end % }}}
	end
end
