%FRICTIONEMULATOR class definition
%
%   Usage:
%      frictionemulator=frictionemulator();

classdef frictionemulator
	properties (SetAccess=public) 
		module_dir = '';
		pt_name = '';
		py_name = '';
		C = NaN;
	end
	methods
		function self = extrude(self,md) % {{{
			self.C=project3d(md,'vector',self.C,'type','node','layer',1);
		end % }}}
		function self = frictionemulator(varargin) % {{{
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{
			md = checkfield(md,'fieldname','friction.module_dir','filepath',1);
			md = checkfield(md,'fieldname','friction.py_name','string',1);
			md = checkfield(md,'fieldname','friction.pt_name','string',1);
			md = checkfield(md,'fieldname','friction.C','timeseries',1,'NaN',1,'Inf',1);
		end % }}}
		function disp(self) % {{{
			disp(sprintf('Basal shear stress parameters for pre-trained python emulator'));
         fielddisplay(self,'module_dir', 'directory of the emulator module');
			fielddisplay(self,'pt_name', 'name of the checkpoint file for pre-trained ML model');
			fielddisplay(self,'py_name', 'name of the python file that defines ML architecture');
			fielddisplay(self,'C','friction coefficient [SI]');
		end % }}}
		function marshall(self,prefix,md,fid) % {{{
			yts=md.constants.yts;
			WriteData(fid,prefix,'name','md.friction.law','data',20,'format','Integer');
			WriteData(fid,prefix,'class','friction','object',self,'fieldname','module_dir','format','String')
			WriteData(fid,prefix,'class','friction','object',self,'fieldname','pt_name','format','String');
			WriteData(fid,prefix,'class','friction','object',self,'fieldname','py_name','format','String');
			WriteData(fid,prefix,'class','friction','object',self,'fieldname','C','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
		end % }}}
		function savemodeljs(self,fid,modelname) % {{{
         error('not implemented yet!');
		end % }}}
	end
end
