%SOLIDEARTH class definition
%
%   Usage:
%      solidearth=solidearth();
%      solidearth=solidearth('earth');

classdef solidearth
	properties (SetAccess=public) 
		settings          = solidearthsettings();
		external          = [];
		lovenumbers       = lovenumbers();
		rotational        = rotational();
		planetradius      = planetradius('earth');
		requested_outputs = {};
		transitions       = {};
		transfercount     = [];
		partitionice      = [];
		partitionhydro    = [];
		partitionocean    = [];
	end
	methods (Static)
		function self = loadobj(self) % {{{
			% This function is directly called by matlab when a model object is
			% loaded. If the input is a struct it is an old version of this class and
			% old fields must be recovered (make sure they are in the deprecated
			% model properties)

			if isstruct(self)
				% 2021, Jan 10
				if isfield(self,'sealevel')
					self.initialsealevel = self.sealevel;
				end
				self = structtoobj(solidearth(),self);

			end
		end% }}}
	end
	methods
		function self = solidearth(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self,'earth');
				case 1
					self=setdefaultparameters(self,varargin{:});
				otherwise
					error('solidearth constructor error message: zero or one argument only!'); 
			end
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   solidearth inputs, forcings and settings:'));

			fielddisplay(self,'planetradius','planet radius [m]');
			fielddisplay(self,'transitions','indices into parts of the mesh that will be icecaps');
			fielddisplay(self,'transfercount','number of icecaps vertices are part of');
			fielddisplay(self,'requested_outputs','additional outputs requested');
			fielddisplay(self,'partitionice','ice partition vector for barystatic contribution');
			fielddisplay(self,'partitionhydro','hydro partition vector for barystatic contribution');
			fielddisplay(self,'partitionocean','ocean partition vector for barystatic contribution');
			if isempty(self.external), fielddisplay(self,'external','external solution, of the type solidearthsolution'); end
			self.settings.disp();
			self.lovenumbers.disp();
			self.rotational.disp();
			if ~isempty(self.external),
				self.external.disp();
			end

		end % }}}
		function self = setdefaultparameters(self,planet) % {{{

			%output default:
			self.requested_outputs={'default'};

			%transitions should be a cell array of vectors:
			self.transitions={};
			self.transfercount=[0];

			%no partitions requested for barystatic contribution:
			self.partitionice=[];
			self.partitionhydro=[];
			self.partitionocean=[];

			%no external solutions by default:
			self.external=[];

			%planet radius
			self.planetradius= planetradius(planet);
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			if ~ismember('SealevelchangeAnalysis',analyses) | (strcmp(solution,'TransientSolution') & md.transient.isslc==0), 
				return; 
			end

			md = checkfield(md,'fieldname','solidearth.requested_outputs','stringrow',1);

			self.settings.checkconsistency(md,solution,analyses);
			self.lovenumbers.checkconsistency(md,solution,analyses);
			self.rotational.checkconsistency(md,solution,analyses);
			if ~isempty(self.external),
				if ~isa(self.external,'solidearthsolution'),
					error('solidearth consistency check: external field should be a solidearthsolution');
				end
				self.external.checkconsistency(md,solution,analyses);
			end

		end % }}}
		function list=defaultoutputs(self,md) % {{{
			list = {'Sealevel'};
		end % }}}
		function marshall(self,prefix,md,fid) % {{{

			WriteData(fid,prefix,'object',self,'fieldname','planetradius','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','transitions','format','MatArray');
			WriteData(fid,prefix,'object',self,'fieldname','transfercount','format','DoubleMat','mattype',1);

			if ~isempty(self.partitionice),
				npartice=max(self.partitionice)+2;
			else
				npartice=0;
			end

			if ~isempty(self.partitionhydro),
				nparthydro=max(self.partitionhydro)+2;
			else
				nparthydro=0;
			end
			if ~isempty(self.partitionocean),
				npartocean=max(self.partitionocean)+2;
			else
				npartocean=0;
			end

			WriteData(fid,prefix,'object',self,'fieldname','partitionice','mattype',1,'format','DoubleMat');
			WriteData(fid,prefix,'data',npartice,'format','Integer','name','md.solidearth.npartice');
			WriteData(fid,prefix,'object',self,'fieldname','partitionhydro','mattype',1,'format','DoubleMat');
			WriteData(fid,prefix,'data',nparthydro,'format','Integer','name','md.solidearth.nparthydro');
			WriteData(fid,prefix,'object',self,'fieldname','partitionocean','mattype',1,'format','DoubleMat');
			WriteData(fid,prefix,'data',npartocean,'format','Integer','name','md.solidearth.npartocean');

			self.settings.marshall(prefix,md,fid);
			self.lovenumbers.marshall(prefix,md,fid);
			self.rotational.marshall(prefix,md,fid);
			if ~isempty(self.external),
				WriteData(fid,prefix,'data',1,'format','Integer','name','md.solidearth.isexternal');
				self.external.marshall(prefix,md,fid);
			else
				WriteData(fid,prefix,'data',0,'format','Integer','name','md.solidearth.isexternal');
			end

			%process requested outputs
			outputs = self.requested_outputs;
			pos = find(ismember(outputs,'default'));
			if ~isempty(pos),
				outputs(pos) = [];                         %remove 'default' from outputs
				outputs      = [outputs defaultoutputs(self,md)]; %add defaults
			end
			WriteData(fid,prefix,'data',outputs,'name','md.solidearth.requested_outputs','format','StringArray');

		end % }}}
		function self = extrude(self,md) % {{{
		end % }}}
		function savemodeljs(self,fid,modelname) % {{{
		
			self.settings.savemodeljs(fid,modelname);
			self.lovenumbers.savemodeljs(fid,modelname);
			self.rotational.savemodeljs(fid,modelname);
			if ~isempty(self.external),
				self.external.savemodeljs(fid,modelname);
			end
			writejscellstring(fid,[modelname '.solidearth.requested_outputs'],self.requested_outputs);
			writejscellarray(fid,[modelname '.solidearth.transitions'],self.transitions);
			writejscellarray(fid,[modelname '.solidearth.partition'],self.partition);
		end % }}}
	end
end
