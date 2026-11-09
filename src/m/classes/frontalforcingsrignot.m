%FRONTALFORCINGSRIGNOT class definition
%
%   Usage:
%      frontalforcingsrignot=frontalforcingsrignot();

classdef frontalforcingsrignot
	properties (SetAccess=public) 
		basin_id             = NaN;
		num_basins           = 0;
		subglacial_discharge = NaN;
		thermalforcing       = NaN;
		isdischargerunoff    = 0;
	end
	methods
		function self = frontalforcingsrignot(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				case 1
					inputstruct=varargin{1};
					list1 = properties('frontalforcingsrignot');
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
		    % nothing for now
		end % }}}
		function self = setdefaultparameters(self) % {{{

			self.basin_id             = NaN;
			self.num_basins           = 0;
			self.subglacial_discharge = NaN;
			self.thermalforcing       = NaN;
			self.isdischargerunoff    = 0;
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{
			%Early return
			if (~strcmp(solution,'TransientSolution') | md.transient.ismovingfront==0), return; end

			md = checkfield(md,'fieldname','frontalforcings.isdischargerunoff','numel',1,'values',[0 1]);
			md = checkfield(md,'fieldname','frontalforcings.num_basins','numel',1,'NaN',1,'Inf',1,'>',0);
			md = checkfield(md,'fieldname','frontalforcings.basin_id','Inf',1,'>=',0,'<=',md.frontalforcings.num_basins,'size',[md.mesh.numberofelements 1]);
			if self.isdischargerunoff
				if md.transient.issmb==0
					error('frontalforcings.isdischargerunoff requires md.transient.issmb=1');
				end
				%SmbRunoffEnum is only provided by some SMB models
				if ~ismember(class(md.smb),{'SMBpddSicopolis','SMBgemb','SMBsemic','SMBgradientscomponents'})
					error(['frontalforcings.isdischargerunoff not supported with smb class ' class(md.smb)]);
				end
			else
				md = checkfield(md,'fieldname','frontalforcings.subglacial_discharge','>=',0,'NaN',1,'Inf',1,'timeseries',1);
			end
			md = checkfield(md,'fieldname','frontalforcings.thermalforcing','NaN',1,'Inf',1,'timeseries',1);

		end % }}}
		function disp(self) % {{{
			disp(sprintf('   Frontalforcings parameters:'));
			fielddisplay(self,'basin_id','basin ID for elements');
			fielddisplay(self,'num_basins', 'number of basins');
			fielddisplay(self,'subglacial_discharge','sum of subglacial discharge for each basin [m^3/d]');
			fielddisplay(self,'thermalforcing','thermal forcing [∘C]');
			fielddisplay(self,'isdischargerunoff','if 1, discharge = basin-integrated SMB runoff over ice-covered area; subglacial_discharge is ignored');
		end % }}}
		function marshall(self,prefix,md,fid) % {{{
			WriteData(fid,prefix,'name','md.frontalforcings.parameterization','data',2,'format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','basin_id','data',self.basin_id-1,'name','md.frontalforcings.basin_id','format','IntMat','mattype',2); %0-indexed
			WriteData(fid,prefix,'object',self,'fieldname','num_basins','format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','isdischargerunoff','format','Boolean');
			if ~self.isdischargerunoff
				WriteData(fid,prefix,'object',self,'fieldname','subglacial_discharge','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			end
			WriteData(fid,prefix,'object',self,'fieldname','thermalforcing','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
		end % }}}
	end
end
