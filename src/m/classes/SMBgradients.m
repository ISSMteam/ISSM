%SMBgradients Class definition
%
%   Usage:
%      SMBgradients=SMBgradients();

classdef SMBgradients
	properties (SetAccess=public)
		href              = NaN;
		smbref            = NaN;
		b_pos             = NaN;
		b_neg             = NaN;
		steps_per_step    = 1;
		averaging         = 0;
		requested_outputs = {};
	end
	methods
		function self = SMBgradients(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = extrude(self,md) % {{{

			%Nothing for now

		end % }}}
		function list = defaultoutputs(self,md) % {{{
			list = {'SmbMassBalance'};
		end % }}}
		function self = initialize(self,md) % {{{

			%Nothing done for now

		end % }}}
		function self = setdefaultparameters(self) % {{{

			%output default:
			self.requested_outputs={'default'};

		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			if ismember('MasstransportAnalysis',analyses),
				md = checkfield(md,'fieldname','smb.href','timeseries',1,'NaN',1,'Inf',1);
				md = checkfield(md,'fieldname','smb.smbref','timeseries',1,'NaN',1,'Inf',1);
				if max(max(abs(md.smb.smbref(1:end-1,:))))<1
					disp('!!! Warning: SMBgradients now expects smbref to be in m/yr ice eq. instead of mm/yr water eq.');
				end
				md = checkfield(md,'fieldname','smb.b_pos','timeseries',1,'NaN',1,'Inf',1);
				md = checkfield(md,'fieldname','smb.b_neg','timeseries',1,'NaN',1,'Inf',1);
			end
			md = checkfield(md,'fieldname','smb.steps_per_step','>=',1,'numel',[1]);
			md = checkfield(md,'fieldname','smb.averaging','numel',[1],'values',[0 1 2]);
			md = checkfield(md,'fieldname','smb.requested_outputs','stringrow',1);
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   surface forcings parameters:'));

			disp(sprintf('\n   SMB gradients parameters:'));
			fielddisplay(self,'href','reference elevation from which deviation is used to calculate SMB adjustment in smb gradients method [m]');
			fielddisplay(self,'smbref','reference smb from which deviation is calculated in smb gradients method [m/yr ice equiv]');
			fielddisplay(self,'b_pos','slope of hs - smb regression line for accumulation regime required if smb gradients is activated');
			fielddisplay(self,'b_neg','slope of hs - smb regression line for ablation regime required if smb gradients is activated');
			fielddisplay(self, 'steps_per_step', 'number of smb steps per time step');
			fielddisplay(self, 'averaging', 'averaging methods from short to long steps');
			disp(sprintf('%51s  0: Arithmetic (default)',' '));
			disp(sprintf('%51s  1: Geometric',' '));
			disp(sprintf('%51s  2: Harmonic',' '));
			fielddisplay(self,'requested_outputs','additional outputs requested');

		end % }}}
		function marshall(self,prefix,md,fid) % {{{

			yts=md.constants.yts;

			WriteData(fid,prefix,'name','md.smb.model','data',6,'format','Integer');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','href','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','smbref','format','DoubleMat','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','b_pos','format','DoubleMat','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','b_neg','format','DoubleMat','mattype',1,'scale',1./yts, ...
				  'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid, prefix, 'object', self, 'fieldname', 'steps_per_step', 'format', 'Integer');
			WriteData(fid, prefix, 'object', self, 'fieldname', 'averaging', 'format', 'Integer');

			%process requested outputs
			outputs = self.requested_outputs;
			pos  = find(ismember(outputs,'default'));
			if ~isempty(pos),
				outputs(pos) = [];                                %remove 'default' from outputs
				outputs      = [outputs defaultoutputs(self,md)]; %add defaults
			end
			WriteData(fid,prefix,'data',outputs,'name','md.smb.requested_outputs','format','StringArray');

		end % }}}
	end
end
