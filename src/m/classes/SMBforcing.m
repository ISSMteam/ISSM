%SMBforcing Class definition
%
%   Usage:
%      SMB=SMBforcing();

classdef SMBforcing
	properties (SetAccess=public)
		mass_balance      = NaN;
		steps_per_step    = 1;
		requested_outputs = {};
		averaging         = 0;
	end
	methods
		function self = SMBforcing(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				case 1
					inputstruct=varargin{1};
					list1 = properties('SMBforcing');
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
		function list = defaultoutputs(self,md) % {{{
			list = {'SmbMassBalance'};
		end % }}}
		function self = extrude(self,md) % {{{
			self.mass_balance=project3d(md,'vector',self.mass_balance,'type','node');
		end % }}}
		function self = initialize(self,md) % {{{

			if isnan(self.mass_balance)
				self.mass_balance=zeros(md.mesh.numberofvertices,1);
				disp('      no smb.mass_balance specified: values set as zero');
			end

		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{
			if (strcmp(solution,'TransientSolution') & md.transient.issmb == 0), return; end
			if ismember('MasstransportAnalysis',analyses),
				md = checkfield(md,'fieldname','smb.mass_balance','timeseries',1,'NaN',1,'Inf',1);
			end
			if ismember('BalancethicknessAnalysis',analyses),
				md = checkfield(md,'fieldname','smb.mass_balance','size',[md.mesh.numberofvertices 1],'NaN',1,'Inf',1);
			end
			md = checkfield(md,'fieldname','smb.steps_per_step','>=',1,'numel',[1]);
			md = checkfield(md,'fieldname','smb.requested_outputs','stringrow',1);
			md = checkfield(md,'fieldname','smb.averaging','numel',[1],'values',[0 1 2]);
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   surface forcings parameters:'));
			fielddisplay(self,'mass_balance','surface mass balance [m/yr ice eq]');
			fielddisplay(self, 'steps_per_step', 'number of smb steps per time step');
			fielddisplay(self,'requested_outputs','additional outputs requested');
			fielddisplay(self,'averaging','averaging methods from short to long steps');
			disp(sprintf('%51s  0: Arithmetic (default)',' '));
			disp(sprintf('%51s  1: Geometric',' '));
			disp(sprintf('%51s  2: Harmonic',' '));
		end % }}}
		function marshall(self,prefix,md,fid) % {{{

			yts=md.constants.yts;

			WriteData(fid,prefix,'name','md.smb.model','data',1,'format','Integer');
			WriteData(fid,prefix,'object',self,'class','smb','fieldname','mass_balance','format','DoubleMat','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1,'yts',yts);
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
		function savemodeljs(self,fid,modelname) % {{{

			writejs1Darray(fid,[modelname '.smb.mass_balance'],self.mass_balance);
			writejscellstring(fid,[modelname '.smb.requested_outputs'],self.requested_outputs);

		end % }}}
		function self = setdefaultparameters(self) % {{{

			%output default:
			self.requested_outputs={'default'};

		end % }}}
	end
end
