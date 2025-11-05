%CALVINCREVASSEDEPTH class definition
%
%   Usage:
%      calvingcrevassedepth=calvingcrevassedepth();

classdef calvingcrevassedepth
	properties (SetAccess=public) 
		crevasse_opening_stress=1;
		crevasse_threshold     =1.;
		water_height = 0.;
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
						if ismember(fieldname,list2),
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
			
			crevasse_threshold      = 1.;
			crevasse_opening_stress = 1;
         self.water_height       = 0.;
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{
			%Early return
			if (~strcmp(solution,'TransientSolution') | md.transient.ismovingfront==0), return; end

			md = checkfield(md,'fieldname','calving.crevasse_opening_stress','numel',[1],'values',[0,1,2]);
         md = checkfield(md,'fieldname','calving.crevasse_threshold','numel',[1],'>',0,'<=',1);
			md = checkfield(md,'fieldname','calving.water_height','NaN',1,'Inf',1,'timeseries',1,'>=',0);
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   Calving Pi parameters:'));
			fielddisplay(self,'crevasse_opening_stress','0: stress only in the ice-flow direction, 1: max principal, 2: buttressing based');
			fielddisplay(self,'crevasse_threshold','ratio of full thickness to calve (e.g. 0.75 is for 75% of the total ice thickness)');
			fielddisplay(self,'water_height','water height in the crevasse [m]');

		end % }}}
		function marshall(self,prefix,md,fid) % {{{
			yts=md.constants.yts;
			WriteData(fid,prefix,'name','md.calving.law','data',6,'format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','crevasse_opening_stress','format','Integer');
         WriteData(fid,prefix,'object',self,'fieldname','crevasse_threshold','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','water_height','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
		end % }}}
	end
end
