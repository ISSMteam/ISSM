%SMBemulator Class definition
%
%   Usage:
%      SMBemulator=SMBemulator();

classdef SMBemulator
	properties (SetAccess=public)

		mass_balance          = NaN;
		elev                  = NaN;
		al                    = NaN;
		st                    = NaN;
		tt                    = NaN;
		swd                   = NaN;
		lwd                   = NaN;
		swu                   = NaN;
		lwu                   = NaN;
		shf                   = NaN;
		lhf                   = NaN;
		requested_outputs     = {};
      module_dir = '';
		pt_name = '';
		py_name = '';
	end
	methods
		function self = SMBemulator(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = extrude(self,md) % {{{
			self.mass_balance=project3d(md,'vector',self.mass_balance,'type','node');
		end % }}}
		function list = defaultoutputs(self,md) % {{{
			list = {'SmbMassBalance'};
		end % }}}
		function self = initialize(self,md) % {{{
			if isnan(self.mass_balance)
				self.mass_balance=zeros(md.mesh.numberofvertices,1);
				disp('      no smb.mass_balance specified: values set as zero');
			end
		end % }}}
		function self = setdefaultparameters(self) % {{{
			self.requested_outputs={'default'};
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{
			if (strcmp(solution,'TransientSolution') & md.transient.issmb == 0), return; end
			if ismember('MasstransportAnalysis',analyses),
			md = checkfield(md,'fieldname','smb.mass_balance','timeseries',1,'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','smb.elev','timeseries',1,'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','smb.al','timeseries',1,'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','smb.st','timeseries',1,'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','smb.tt','timeseries',1,'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','smb.swd','timeseries',1,'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','smb.lwd','timeseries',1,'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','smb.swu','timeseries',1,'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','smb.lwu','timeseries',1,'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','smb.shf','timeseries',1,'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','smb.lhf','timeseries',1,'NaN',1,'Inf',1);
			end
			md = checkfield(md,'fieldname','smb.requested_outputs','stringrow',1);
			md = checkfield(md,'fieldname','smb.module_dir','filepath',1);
			md = checkfield(md,'fieldname','smb.py_name','string',1);
			md = checkfield(md,'fieldname','smb.pt_name','string',1);
		end % }}}
		function disp(self) % {{{
			disp(sprintf('\n  smb emulator based on MAR-IA :'));
         fielddisplay(self, 'mass_balance', 'surface mass balance for validation purpose [m/yr ice eq]'); 
         fielddisplay(self, 'elev', 'surface elevation for validation purpose');
         fielddisplay(self, 'al', 'albedo');
         fielddisplay(self, 'st', 'surface temperature');
         fielddisplay(self, 'tt', 'two meter air temperature');
         fielddisplay(self, 'swd', 'short wave radiation down');
         fielddisplay(self, 'lwd', 'long wave radiation down');
         fielddisplay(self, 'swu', 'short wave radiation up');
         fielddisplay(self, 'lwu', 'long wave radiation up');
         fielddisplay(self, 'shf', 'sensible heat flux');
         fielddisplay(self, 'lhf', 'latent heat flux');
			fielddisplay(self,'requested_outputs','additional outputs requested');
         fielddisplay(self,'module_dir', 'directory of the emulator module');
			fielddisplay(self,'pt_name', 'name of the checkpoint file for pre-trained ML model');
			fielddisplay(self,'py_name', 'name of the python file that defines ML architecture');
		end % }}}
		function marshall(self,prefix,md,fid) % {{{
			WriteData(fid,prefix,'name','md.smb.model','data',20,'format','Integer');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','mass_balance','format','DoubleMat','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1,'yts',yts); % double check units
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','elev','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',yts); % unit of elev?
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','al','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',yts);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','st','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',yts);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','tt','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',yts);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','swd','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',yts);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','lwd','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',yts);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','swu','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',yts);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','lwu','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',yts);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','shf','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',yts);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','lhf','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',yts);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','module_dir','format','String')
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','pt_name','format','String');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','py_name','format','String');

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
