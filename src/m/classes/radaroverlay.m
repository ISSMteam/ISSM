%RADAROVERLAY class definition
%
%   Usage:
%      radaroverlay=radaroverlay();

classdef radaroverlay
	properties (SetAccess=public) 
		pwr = NaN;
		x   = NaN;
		y   = NaN;
		outerindex = NaN;
		outerx = NaN;
		outery = NaN;
		outerheight = NaN;
	end
	methods
		function self = radaroverlay(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = setdefaultparameters(self) % {{{

		end % }}}
		function disp(self) % {{{
			disp(sprintf('   radaroverlay parameters:'));

			fielddisplay(self,'pwr','radar power image (matrix)');
			fielddisplay(self,'x','corresponding x coordinates [m]');
			fielddisplay(self,'y','corresponding y coordinates [m]');
			fielddisplay(self,'outerindex','outer triangulation corresponding to space between mesh and the bounding box');
			fielddisplay(self,'outerx','outer x corresponding to space between mesh and the bounding box');
			fielddisplay(self,'outery','outer y corresponding to space between mesh and the bounding box');
			fielddisplay(self,'outerheight','outer height corresponding to space between mesh and the bounding box');

		end % }}}
		function savemodeljs(self,fid,modelname) % {{{
		
			if ~isnan(self.pwr),
				writejs1Darray(fid,[modelname '.radaroverlay.xlim'],[min(self.x) max(self.x)]);
				writejs1Darray(fid,[modelname '.radaroverlay.ylim'],[min(self.y) max(self.y)]);
				writejs2Darray(fid,[modelname '.radaroverlay.outerindex'],self.outerindex);
				writejs1Darray(fid,[modelname '.radaroverlay.outerx'],self.outerx);
				writejs1Darray(fid,[modelname '.radaroverlay.outery'],self.outery);
				writejs1Darray(fid,[modelname '.radaroverlay.outerheight'],self.outerheight)
			end

		end % }}}
	end
end
