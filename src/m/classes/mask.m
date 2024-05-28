%MASK class definition
%
%   Usage:
%      mask=mask();

classdef mask
	properties (SetAccess=public) 
		ocean_levelset = NaN;
		ice_levelset   = NaN;
	end
	methods (Static)
		function self = loadobj(self) % {{{
			% This function is directly called by matlab when a model object is
			% loaded. Update old properties here

			%2014 February 5th. NOT NEED ANYMORE, AND CAN LEAD TO CONFUSION for sea-level models. 
			%if numel(self.ice_levelset)>1 & all(self.ice_levelset>=0),
			%	disp('WARNING: md.mask.ice_levelset>=0, you probably need to change the sign of this levelset');
			%end

			%2020 May 15th
			if isstruct(self)
				selfnew = mask();
				if isfield(self,'ice_levelset')
					selfnew.ice_levelset = self.ice_levelset;
				end
				if isfield(self,'ocean_levelset')
					selfnew.ocean_levelset = self.ocean_levelset;
				end
				if isfield(self,'groundedice_levelset')
					disp('WARNING: automatically updated md.mask as groundedice_levelset is now ocean_levelset');
					selfnew.ocean_levelset = self.groundedice_levelset;
				end
				self = selfnew;
			end
		end % }}}
	end
	methods
		function disp(self) % {{{
			disp(sprintf('   masks:'));

			fielddisplay(self,'ocean_levelset','presence of ocean if < 0, coastline/grounding line if = 0, no ocean if > 0');
			fielddisplay(self,'ice_levelset','presence of ice if < 0, icefront position if = 0, no ice if > 0');
		end % }}}
		function self = setdefaultparameters(self) % {{{
		end % }}}
		function self = extrude(self,md) % {{{
			self.ocean_levelset=project3d(md,'vector',self.ocean_levelset,'type','node');
			self.ice_levelset=project3d(md,'vector',self.ice_levelset,'type','node');
		end % }}}
		function self = mask(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self=setocean(self,varargin) % {{{
			if nargin==3,
				ocean=varargin{1}; 
				index=varargin{2};
				self.ocean_levelset(index)=-ocean;
			elseif nargin==2,
				ocean=varargin{1}; 
				self.ocean_levelset=-ocean;
			else error('oceanset error message: not supported yet');
			end
		end % }}}
		function self=setice(self,varargin) % {{{
			if nargin==3,
				ice=varargin{1}; 
				index=varargin{2};
				self.ice_levelset(index)=-ice;
			elseif nargin==2,
				ice=varargin{1}; 
				self.ice_levelset=-ice;
			else error('iceset error message: not supported yet');
			end

		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{
			if strcmpi(solution,'LoveSolution'), return; end;

			if isa(self.ocean_levelset,'cell'),
				for i=1:length(self.ocean_levelset),
					md = checkfield(md,'field',self.ocean_levelset{i},'NaN',0,'Inf',1,'timeserieslength',1,'Inf',1);
				end
			else
				md = checkfield(md,'fieldname','mask.ocean_levelset','timeseries',1,'NaN',1);
			end
			
			if isa(self.ice_levelset,'cell'),
				for i=1:length(self.ice_levelset),
					md = checkfield(md,'field',self.ice_levelset{i},'NaN',0,'Inf',1,'timeserieslength',1,'Inf',1);
				end
			else
				md = checkfield(md,'fieldname','mask.ice_levelset','timeseries',1,'NaN',1);
				isice=(md.mask.ice_levelset<=0);
				if sum(isice)==0,
					warning('no ice present in the domain');
				end
			end
		end % }}}
		function marshall(self,prefix,md,fid) % {{{

			if isa(self.ocean_levelset,'cell'),
				WriteData(fid,prefix,'object',self,'fieldname','ocean_levelset','name','md.mask.ocean_levelset','format','MatArray','timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			else
				WriteData(fid,prefix,'object',self,'fieldname','ocean_levelset','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			end

			if isa(self.ice_levelset,'cell'),
				WriteData(fid,prefix,'object',self,'fieldname','ice_levelset','name','md.mask.ice_levelset','format','MatArray','timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			else
				WriteData(fid,prefix,'object',self,'fieldname','ice_levelset','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			end
		end % }}}
		function savemodeljs(self,fid,modelname) % {{{
		
			writejs1Darray(fid,[modelname '.mask.ocean_levelset'],self.ocean_levelset);
			writejs1Darray(fid,[modelname '.mask.ice_levelset'],self.ice_levelset);

		end % }}}
	end
end
