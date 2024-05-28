%TOOLKITS class definition
%
%   Usage:
%      self=toolkits();

classdef toolkits < dynamicprops
	properties (SetAccess=public) 
		DefaultAnalysis  = struct();
		RecoveryAnalysis = struct();
		%The other properties are dynamic
	end
	methods (Static)
		function self = loadobj(self) % {{{
			% This function is directly called by matlab when a model object is
			% loaded. Update old properties here

			if isempty(fieldnames(self.RecoveryAnalysis));
				disp('WARNING: updating toolkits (RecoveryAnalysis not set)');
				self.RecoveryAnalysis  = self.DefaultAnalysis;
			end
		end% }}}
	end
	methods
		function self = toolkits(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				case 1
					self=structtoobj(self,varargin{1});
				otherwise
					error('constructor not supported');
				end
			end % }}}
		function self = addoptions(self,analysis,varargin) % {{{
		%ADDOPTIONS - add analysis to md.toolkits.analysis
		%
		%   Optional third parameter adds toolkits options to analysis.
		%
		%   Usage:
		%      md.toolkits=addoptions(md.toolkits,'StressbalanceAnalysis',FSoptions());
		%      md.toolkits=addoptions(md.toolkits,'StressbalanceAnalysis');

			%Create dynamic property if property does not exist yet
			if ~ismember(analysis,properties(self)),
				self.addprop(analysis);
			end

			%Add toolkits options to analysis
			if nargin==3,
				self.(analysis) = varargin{1};
			end
		end
		%}}}
		function self = setdefaultparameters(self) % {{{

			%default toolkits: 
			if IssmConfig('_HAVE_PETSC_'),
				%MUMPS is the default toolkits
				if IssmConfig('_HAVE_MUMPS_'),
					self.DefaultAnalysis           = mumpsoptions();
				else
					self.DefaultAnalysis           = iluasmoptions(); 
				end
			else
				if IssmConfig('_HAVE_MUMPS_'),
					self.DefaultAnalysis           = issmmumpssolver(); 
				elseif IssmConfig('_HAVE_GSL_'),
					self.DefaultAnalysis           = issmgslsolver(); 
				else 
					disp('WARNING: Need at least MUMPS or GSL to define an ISSM solver type, no default solver assigned');
				end
			end

			%Use same solver for Recovery mode
			self.RecoveryAnalysis = self.DefaultAnalysis;


		end % }}}
		function disp(self) % {{{
			analyses=properties(self);
			disp(sprintf('List of toolkits options per analysis:\n'));
			for i=1:numel(analyses),
				analysis=analyses{i};
				disp([analysis ':']);
				disp(self.(analysis));
			end
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{
			analyses=properties(self);
			for i=1:numel(analyses),
				switch analyses{i}
					case 'DefaultAnalysis'
					case 'RecoveryAnalysis'
					case 'StressbalanceAnalysis'
					case 'StressbalanceVerticalAnalysis'
					case 'GLheightadvectionAnalysis'
					case 'MasstransportAnalysis'
					case 'ThermalAnalysis'
					case 'EnthalpyAnalysis'
					case 'AdjointBalancethicknessAnalysis'
					case 'BalancethicknessAnalysis'
					case 'Balancethickness2Analysis'
					case 'BalancethicknessSoftAnalysis'
					case 'BalancevelocityAnalysis'
					case 'DamageEvolutionAnalysis'
					case 'LoveAnalysis'
					case 'EsaAnalysis'
					case 'SealevelchangeAnalysis'
					case 'FreeSurfaceBaseAnalysis'
					case 'FreeSurfaceTopAnalysis'
					case 'LevelsetAnalysis'
					case 'DebrisAnalysis'
					case 'L2ProjectionBaseAnalysis'
					case 'ExtrudeFromBaseAnalysis'
					case 'ExtrudeFromTopAnalysis'
					otherwise
						md = checkmessage(md,['md.toolkits.' analyses{i} ' not supported yet']);
				end
				if isempty(fieldnames(self.(analyses{i})))
					md = checkmessage(md,['md.toolkits.' analyses{i} ' is empty']);
				end
			end
		end % }}}
		function ToolkitsFile(toolkits,filename) % {{{
			%TOOLKITSFILE - build toolkits file
			%
			%   Build a Petsc compatible options file, from the toolkits model field and return options string.
			%   This file will also be used when the toolkit used is 'issm' instead of 'petsc'.
			%
			%   Usage:     ToolkitsFile(toolkits,filename);

			%open file for writing
			fid=fopen(filename,'w');
			if fid==-1,
				error(['ToolkitsFile error: could not open ' filename ' for writing']);
			end

			%write header
			fprintf(fid,'%s%s%s\n','%Toolkits options file: ',filename,' written from Matlab toolkits array');

			%start writing options
			analyses=properties(toolkits);
			for i=1:numel(analyses),
				analysis=analyses{i};
				options=toolkits.(analysis);

				%first write analysis:
				fprintf(fid,'\n+%s\n',analysis); %append a + to recognize it's an analysis enum

				%now, write options
				optionslist=fieldnames(options);
				for j=1:numel(optionslist),
					optionname=optionslist{j};
					optionvalue=options.(optionname);

					if isempty(optionvalue),
						%this option has only one argument
						fprintf(fid,'-%s\n',optionname);
					else
						%option with value. value can be string or scalar
						if isnumeric(optionvalue),
							fprintf(fid,'-%s %g\n',optionname,optionvalue);
						elseif ischar(optionvalue),
							fprintf(fid,'-%s %s\n',optionname,optionvalue);
						else
							error(['ToolkitsFile error: option ' optionname ' is not well formatted']);
						end
					end
				end
			end

			fclose(fid);
		end %}}}
		function savemodeljs(self,fid,modelname) % {{{

			analyses=properties(self);
			for i=1:numel(analyses),
				if isempty(fieldnames(self.(analyses{i})))
					error(['md.toolkits.' analyses{i} ' is empty']);
				else
					writejsstruct(fid,[modelname '.toolkits.' analyses{i}],self.(analyses{i}));
				end
			end

		end % }}}
	end
 end
