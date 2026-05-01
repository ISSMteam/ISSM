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
		steps_per_step        = 1;
		requested_outputs     = {};
		averaging             = 0;
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
			self.al=project3d(md,'vector',self.al,'type','node');
			self.st=project3d(md,'vector',self.st,'type','node');
			self.tt=project3d(md,'vector',self.tt,'type','node');
			self.swd=project3d(md,'vector',self.swd,'type','node');
			self.lwd=project3d(md,'vector',self.lwd,'type','node');
			self.swu=project3d(md,'vector',self.swu,'type','node');
			self.lwu=project3d(md,'vector',self.lwu,'type','node');
			self.shf=project3d(md,'vector',self.shf,'type','node');
			self.lhf=project3d(md,'vector',self.lhf,'type','node');
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
			md = checkfield(md,'fieldname','smb.steps_per_step','>=',1,'numel',[1]);
			md = checkfield(md,'fieldname','smb.requested_outputs','stringrow',1);
			md = checkfield(md,'fieldname','smb.averaging','numel',[1],'values',[0 1 2]);
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
			fielddisplay(self, 'steps_per_step', 'number of smb steps per time step');
			fielddisplay(self,'requested_outputs','additional outputs requested');
			fielddisplay(self,'averaging','averaging methods from short to long steps');
         fielddisplay(self,'module_dir', 'directory of the emulator module');
			fielddisplay(self,'pt_name', 'name of the checkpoint file for pre-trained ML model');
			fielddisplay(self,'py_name', 'name of the python file that defines ML architecture');
		end % }}}
		function marshall(self,prefix,md,fid) % {{{
			yts = md.constants.yts;
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
			WriteData(fid,prefix,'object',self,'fieldname','steps_per_step','format','Integer');
			 WriteData(fid,prefix,'object',self,'fieldname','averaging','format','Integer');

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
