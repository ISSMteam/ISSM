%RIFTS class definition
%
%   Usage:
%      rifts=rifts();

classdef rifts
	properties (SetAccess=public) 
		riftstruct     = NaN;
		riftproperties = NaN;
	end
	methods
		function self = rifts(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = setdefaultparameters(self) % {{{

		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{
			if isempty(self.riftstruct) | isnans(self.riftstruct),
				numrifts=0;
			else
				numrifts=numel(self.riftstruct);
			end
			if numrifts,
				if ~(strcmp(domaintype(md.mesh),'2Dhorizontal')),
					md = checkmessage(md,['models with rifts are only supported in 2d for now!']);
				end
				if ~isstruct(self.riftstruct),
					md = checkmessage(md,['rifts.riftstruct should be a structure!']);
				end
				if ~isempty(find(md.mesh.segmentmarkers>=2)),
					%We have segments with rift markers, but no rift structure!
					md = checkmessage(md,['model should be processed for rifts (run meshprocessrifts)!']);
				end
				for i=1:numrifts,
					md = checkfield(md,'fieldname',sprintf('rifts.riftstruct(%d).fill',i),'values',{'Air','Ice','Melange','Water'});
				end
			else
				if ~isnans(self.riftstruct),
					md = checkmessage(md,['riftstruct should be NaN since numrifts is 0!']);
				end
			end
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   rifts parameters:'));

			fielddisplay(self,'riftstruct','structure containing all rift information (vertices coordinates, segments, type of melange, ...)');
			fielddisplay(self,'riftproperties','');

		end % }}}
		function marshall(self,prefix,md,fid) % {{{

			%Process rift info
			if isempty(self.riftstruct) | isnans(self.riftstruct),
				numrifts=0;
			else
				numrifts=numel(self.riftstruct);
			end

			numpairs=0;
			for i=1:numrifts,
				numpairs=numpairs+size(self.riftstruct(i).penaltypairs,1);
			end
			
			for i=1:numrifts
				if (strcmpi(self.riftstruct(i).fill,'Air'))
					self.riftstruct(i).fill = 0;
				elseif (strcmpi(self.riftstruct(i).fill,'Ice'))
					self.riftstruct(i).fill = 1;
				elseif (strcmpi(self.riftstruct(i).fill,'Melange'))
					self.riftstruct(i).fill = 2;
				elseif (strcmpi(self.riftstruct(i).fill,'Water'))
					self.riftstruct(i).fill = 3;
				else
					error(['Could not convert string in riftstruct to integer for marshalling']);	
				end
			end

			% 2 for nodes + 2 for elements+ 2 for  normals + 1 for length + 1 for fill + 1 for friction + 1 for fraction + 1 for fractionincrement + 1 for state.
			data=zeros(numpairs,12);
			count=1;
			for i=1:numrifts,
				numpairsforthisrift=size(self.riftstruct(i).penaltypairs,1);
				data(count:count+numpairsforthisrift-1,1:7)=self.riftstruct(i).penaltypairs;
				data(count:count+numpairsforthisrift-1,8)=self.riftstruct(i).fill;
				data(count:count+numpairsforthisrift-1,9)=self.riftstruct(i).friction;
				data(count:count+numpairsforthisrift-1,10)=self.riftstruct(i).fraction;
				data(count:count+numpairsforthisrift-1,11)=self.riftstruct(i).fractionincrement;
				data(count:count+numpairsforthisrift-1,12)=self.riftstruct(i).state;
				count=count+numpairsforthisrift;
			end

			WriteData(fid,prefix,'data',numrifts,'name','md.rifts.numrifts','format','Integer');
			WriteData(fid,prefix,'data',data,    'name','md.rifts.riftstruct','format','DoubleMat','mattype',3);
		end % }}}
		function savemodeljs(self,fid,modelname) % {{{
	
			if isempty(self.riftstruct) | isnans(self.riftstruct),
				numrifts=0;
			else
				numrifts=numel(self.riftstruct);
			end
			
			if numrifts,
				error('rifts savemodeljs error message: not supported yet!');
			end
	
		end % }}}
	end
end
