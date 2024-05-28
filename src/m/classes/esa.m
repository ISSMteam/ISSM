%ESA class definition
%
%   Usage:
%      esa=esa();

classdef esa
	properties (SetAccess=public) 
		deltathickness = NaN;
		love_h         = 0; %provided by PREM model
		love_l         = 0; %ideam
		degacc         = 0;
		hemisphere		= 0;
		requested_outputs      = {};
		transitions    = {};
	end
	methods
		function self = esa(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = setdefaultparameters(self) % {{{
		
		%numerical discretization accuracy
		self.degacc=.01;
	
		%computational flags:
		self.hemisphere=0;

		%output default:
		self.requested_outputs={'default'};

		%transitions should be a cell array of vectors: 
		self.transitions={};
		
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			if ~ismember('EsaAnalysis',analyses), return; end
			md = checkfield(md,'fieldname','esa.deltathickness','NaN',1,'Inf',1,'size',[md.mesh.numberofelements 1]);
			md = checkfield(md,'fieldname','esa.love_h','NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','esa.love_l','NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','esa.hemisphere','NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','esa.degacc','size',[1 1],'>=',1e-10);
			md = checkfield(md,'fieldname','esa.requested_outputs','stringrow',1);

			%check that love numbers are provided at the same level of accuracy: 
			if (size(self.love_h,1)~=size(self.love_l,1)),
				error('esa error message: love numbers should be provided at the same level of accuracy');
			end

			%cross check that whereever we have an ice load, the mask is <0 on each vertex: 
			pos=find(self.deltathickness);
			maskpos=md.mask.ice_levelset(md.mesh.elements(pos,:)); 
			[els,vertices]=find(maskpos>0);
			if length(els),
				error('esa checkconsistency fail: there are elements with ice loads where some vertices are not on the ice!');
			end

		end % }}}
		function list=defaultoutputs(self,md) % {{{
			list = {'EsaUmotion'};
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   esa parameters:'));

			fielddisplay(self,'deltathickness','thickness change: ice height equivalent [m]');
			fielddisplay(self,'love_h','load Love number for radial displacement');
			fielddisplay(self,'love_l','load Love number for horizontal displacements');
			fielddisplay(self,'hemisphere','North-south, East-west components of 2-D horiz displacement vector: -1 south, 1 north'); 
			fielddisplay(self,'degacc','accuracy (default .01 deg) for numerical discretization of the Green''s functions');
			fielddisplay(self,'transitions','indices into parts of the mesh that will be icecaps');
			fielddisplay(self,'requested_outputs','additional outputs requested (e.g., EsaUmotion, EsaStrainratexx, EsaRotationrate)');

		end % }}}
		function marshall(self,prefix,md,fid) % {{{
			WriteData(fid,prefix,'object',self,'fieldname','deltathickness','format','DoubleMat','mattype',2);
			WriteData(fid,prefix,'object',self,'fieldname','love_h','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'fieldname','love_l','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'fieldname','hemisphere','format','Integer'); 
			WriteData(fid,prefix,'object',self,'fieldname','degacc','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','transitions','format','MatArray');
			
			%process requested outputs
			outputs = self.requested_outputs;
			pos  = find(ismember(outputs,'default'));
			if ~isempty(pos),
				outputs(pos) = [];                         %remove 'default' from outputs
				outputs      = [outputs defaultoutputs(self,md)]; %add defaults
			end
			WriteData(fid,prefix,'data',outputs,'name','md.esa.requested_outputs','format','StringArray');

		end % }}}
		function savemodeljs(self,fid,modelname) % {{{
		
			writejs1Darray(fid,[modelname '.esa.deltathickness'],self.deltathickness);
			writejs1Darray(fid,[modelname '.esa.love_h'],self.love_h);
			writejs1Darray(fid,[modelname '.esa.love_l'],self.love_l);
			writejsdouble(fid,[modelname '.esa.hemisphere'],self.hemisphere); 
			writejsdouble(fid,[modelname '.esa.degacc'],self.degacc);
			writejscellstring(fid,[modelname '.esa.requested_outputs'],self.requested_outputs);
			writejscellarray(fid,[modelname '.esa.transitions'],self.transitions);
		end % }}}
	end
end
