%SMBhenning Class definition
%
%   Usage:
%      SMBhenning=SMBhenning();

classdef SMBhenning
	properties (SetAccess=public)
		smbref = NaN;
		steps_per_step=1;
		averaging=0;
		requested_outputs      = {};
	end
	methods
		function self = SMBhenning(varargin) % {{{
			switch nargin
				case 0
				case 1
					inputstruct=varargin{1};
					list1 = properties('SMBhenning');
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

			self.smbref=project3d(md,'vector',self.smbref,'type','node');

		end % }}}
		function list = defaultoutputs(self,md) % {{{
			list = {'SmbMassBalance'};
		end % }}}
		function self = initialize(self,md) % {{{

			if isnan(self.smbref)
				self.smbref=zeros(md.mesh.numberofvertices,1);
				disp('      no smb.smbref specified: values set as zero');
			end

		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			if ismember('MasstransportAnalysis',analyses),
				md = checkfield(md,'fieldname','smb.smbref','timeseries',1,'NaN',1,'Inf',1);
			end
			if ismember('BalancethicknessAnalysis',analyses),
				md = checkfield(md,'fieldname','smb.smbref','size',[md.mesh.numberofvertices 1],'NaN',1,'Inf',1);
			end
			md = checkfield(md,'fieldname','smb.steps_per_step','>=',1,'numel',[1]);
			md = checkfield(md,'fieldname','smb.averaging','numel',[1],'values',[0 1 2]);
			md = checkfield(md,'fieldname','smb.requested_outputs','stringrow',1);
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   surface forcings parameters:'));
			fielddisplay(self,'smbref','reference smb from which deviation is calculated [m/yr ice eq]');
			fielddisplay(self, 'steps_per_step', 'number of smb steps per time step');
			fielddisplay(self, 'averaging', 'averaging methods from short to long steps');
			disp(sprintf('%51s  0: Arithmetic (default)',' '));
			disp(sprintf('%51s  1: Geometric',' '));
			disp(sprintf('%51s  2: Harmonic',' '));
			fielddisplay(self,'requested_outputs','additional outputs requested');
		end % }}}
		function marshall(self,prefix,md,fid) % {{{

			yts=md.constants.yts;

			WriteData(fid,prefix,'name','md.smb.model','data',7,'format','Integer');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','smbref','format','DoubleMat','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
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
