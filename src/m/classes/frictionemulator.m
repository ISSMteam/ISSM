%FRICTIONEMULATOR class definition
%
%   Usage:
%      frictionemulator=frictionemulator();

classdef frictionemulator
	properties (SetAccess=public) 
		pt_path = '';
	end
	methods
		function self = extrude(self,md) % {{{
		end % }}}
		function self = frictionemulator(varargin) % {{{
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			md = checkfield(md,'fieldname','friction.pt_path','filepath',1);
		end % }}}
		function disp(self) % {{{
			disp(sprintf('Basal shear stress parameters for pre-trained python emulator'));
			fielddisplay(self,'pt_path', 'path to checkpoint file for pre-trained ML model');
		end % }}}
		function marshall(self,prefix,md,fid) % {{{
			yts=md.constants.yts;

			WriteData(fid,prefix,'name','md.friction.law','data',20,'format','Integer');
			WriteData(fid,prefix,'class','friction','object',self,'fieldname','pt_path','format','String');
		end % }}}
		function savemodeljs(self,fid,modelname) % {{{
         error('not implemented yet!');
		end % }}}
	end
end
