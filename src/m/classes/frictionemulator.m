%FRICTIONEMULATOR class definition
%
%   Usage:
%      frictionemulator=frictionemulator();

classdef frictionemulator
	properties (SetAccess=public) 
		module_dir = '';
		pt_name = '';
		py_name = '';
	end
	methods
		function self = extrude(self,md) % {{{
		end % }}}
		function self = frictionemulator(varargin) % {{{
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{
			md = checkfield(md,'fieldname','friction.module_dir','filepath',1);
			md = checkfield(md,'fieldname','friction.py_name','stringrow',1);
			md = checkfield(md,'fieldname','friction.pt_name','stringrow',1);
		end % }}}
		function disp(self) % {{{
			disp(sprintf('Basal shear stress parameters for pre-trained python emulator'));
         fielddisplay(self,'module_dir', 'directory of the emulator module');
			fielddisplay(self,'pt_name', 'name of the checkpoint file for pre-trained ML model');
			fielddisplay(self,'py_name', 'name of the python file that defines ML architecture');
		end % }}}
		function marshall(self,prefix,md,fid) % {{{
			yts=md.constants.yts;
			WriteData(fid,prefix,'name','md.friction.law','data',20,'format','Integer');
			WriteData(fid,prefix,'class','friction','object',self,'fieldname','module_dir','format','String')
			WriteData(fid,prefix,'class','friction','object',self,'fieldname','pt_name','format','String');
			WriteData(fid,prefix,'class','friction','object',self,'fieldname','py_name','format','String');
		end % }}}
		function savemodeljs(self,fid,modelname) % {{{
         error('not implemented yet!');
		end % }}}
	end
end
