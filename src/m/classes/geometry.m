%GEOMETRY class definition
%
%   Usage:
%      geometry=geometry();

classdef geometry
	properties (SetAccess=public) 
		surface           = NaN;
		thickness         = NaN;
		base              = NaN;
		bed               = NaN;
		hydrostatic_ratio = NaN;
	end
	methods (Static)
		function self = loadobj(self) % {{{
			% This function is directly called by matlab when a model object is
			% loaded. Update old properties here

			%2014 March 26th
			if isstruct(self),
				disp('WARNING: updating geometry');
				disp('         md.geometry.bed        is now md.geometry.base');
				disp('         md.geometry.bathymetry is now md.geometry.bed');
				obj2 = self;
				self = geometry();
				self.surface    = obj2.surface;
				self.thickness  = obj2.thickness;
				self.base       = obj2.bed;
				self.bed        = obj2.bathymetry;
			end

		end% }}}
	end
	methods
		function self = geometry(varargin) % {{{
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

			if strcmpi(solution,'LoveSolution'),
				return; 
			else
				md = checkfield(md,'fieldname','geometry.surface' ,'NaN',1,'Inf',1,'size',[md.mesh.numberofvertices 1]);
				md = checkfield(md,'fieldname','geometry.base'      ,'NaN',1,'Inf',1,'size',[md.mesh.numberofvertices 1]);
				md = checkfield(md,'fieldname','geometry.thickness','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices 1],'>=',0);
				if any(abs(self.thickness-self.surface+self.base)>10^-9),
					md = checkmessage(md,['equality thickness=surface-base violated']);
				end 
				if strcmp(solution,'TransientSolution') & md.transient.isgroundingline,
					md = checkfield(md,'fieldname','geometry.bed','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices 1]);
					if any(self.bed-self.base>10^-12),
						md = checkmessage(md,['base<bed on one or more vertices']);
					end 
					pos = find(md.mask.ocean_levelset>0);
					if any(abs(self.bed(pos)-self.base(pos))>10^-9),
						md = checkmessage(md,['equality base=bed on grounded ice violated']);
					end 
					md = checkfield(md,'fieldname','geometry.bed','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices 1]);
				end
			end
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   geometry parameters:'));

			fielddisplay(self,'surface','ice upper surface elevation [m]');
			fielddisplay(self,'thickness','ice thickness [m]');
			fielddisplay(self,'base','ice base elevation [m]');
			fielddisplay(self,'bed','bed elevation [m]');

		end % }}}
		function marshall(self,prefix,md,fid) % {{{
			length_thickness=size(self.thickness,1);
			if length_thickness==md.mesh.numberofvertices | length_thickness==md.mesh.numberofvertices+1,
				WriteData(fid,prefix,'object',self,'fieldname','thickness','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			elseif length_thickness==md.mesh.numberofelements | length_thickness==md.mesh.numberofelements+1,
				WriteData(fid,prefix,'object',self,'fieldname','thickness','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofelements+1,'yts',md.constants.yts);
			else
				error('geometry thickness time series should be a vertex or element time series');
			end

			WriteData(fid,prefix,'object',self,'fieldname','surface','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'fieldname','base','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'fieldname','bed','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'fieldname','hydrostatic_ratio','format','DoubleMat','mattype',1);
		end % }}}
		function self = extrude(self,md) % {{{
			self.surface=project3d(md,'vector',self.surface,'type','node');
			self.thickness=project3d(md,'vector',self.thickness,'type','node');
			self.hydrostatic_ratio=project3d(md,'vector',self.hydrostatic_ratio,'type','node');
			self.base=project3d(md,'vector',self.base,'type','node');
			self.bed=project3d(md,'vector',self.bed,'type','node');
		end % }}}
		function savemodeljs(self,fid,modelname) % {{{
		
			writejs1Darray(fid,[modelname '.geometry.surface'],self.surface);
			writejs1Darray(fid,[modelname '.geometry.thickness'],self.thickness);
			writejs1Darray(fid,[modelname '.geometry.base'],self.base);
			writejs1Darray(fid,[modelname '.geometry.bed'],self.bed);
			writejs1Darray(fid,[modelname '.geometry.hydrostatic_ratio'],self.hydrostatic_ratio);

		end % }}}
	end
end
