%AMRNEOPZ Class definition
%
%   Usage:
%      md.amr=amrneopz();

classdef amrneopz
	properties (SetAccess=public) 
		level_max								= 0; 
		gradation								= 0;
      lag										= 0;
		groundingline_distance				= 0;
      icefront_distance						= 0;
      thicknesserror_threshold			= 0;
      thicknesserror_groupthreshold 	= 0;
      thicknesserror_maximum				= 0;
		deviatoricerror_threshold			= 0;
		deviatoricerror_groupthreshold	= 0;
		deviatoricerror_maximum				= 0;
		restart									= 0;
	end
   methods (Static)
      %function self = loadobj(self) % {{{
         % This function is directly called by matlab when a model object is
         % loaded. Update old properties here

			%if verLessThan('matlab','7.9'),
         %   disp('Warning: your matlab version is old and there is a risk that load does not work correctly');
         %   disp('         if the model is not loaded correctly, rename temporarily loadobj so that matlab does not use it');

            % This is a Matlab bug: all the fields of md have their default value
            % Example of error message:
            % Warning: Error loading an object of class 'model':
            % Undefined function or method 'exist' for input arguments of type 'cell'
            %
            % This has been fixed in MATLAB 7.9 (R2009b) and later versions
			%end

         %2017 September 15th
         %if isstruct(self),
         %   disp('WARNING: updating amr');
         %   disp('         md.amr.region_level_max is now md.amr.radius_level_max');
         %   disp('         md.amr.region_level_1 is not being used; now gradation is used instead');
         %   obj2 						= self;
         %   self 						= amr();
         %   %Converting region_level_1 to gradation
			%	if(obj2.region_level_max>0 && obj2.region_level_1>obj2.region_level_max) 
			%		alpha=0;
			%		if(obj2.level_max>1) alpha=log(obj2.region_level_1/obj2.region_level_max)/(obj2.level_max-1);end
			%		self.radius_level_max			= obj2.region_level_max;
			%		self.level_max						= obj2.level_max;
			%		self.gradation						= exp(alpha);
         %   	self.lag								= 1.0;
			%		self.groundingline_distance 	= 0;
			%		self.icefront_distance 			= 0;
			%		self.thicknesserror_threshold = 0;
			%		self.deviatoricerror_threshold= 0;
			%	end
         %end

			%2017 November 24th
			%radius_level_max was deleted!

			%2017 December 18th
			%group threshold was inserted

      %end% }}}
   end
	methods
		function self = amrneopz(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = setdefaultparameters(self) % {{{

			%filter parameters:
			self.level_max								= 2;
			self.gradation								= 1.5;
			self.lag										= 1.1;
		
 			%other criterias
         self.groundingline_distance			= 10000;
         self.icefront_distance					= 0;
         self.thicknesserror_threshold			= 0;
         self.thicknesserror_groupthreshold	= 0;
			self.thicknesserror_maximum			= 0; 
			self.deviatoricerror_threshold		= 0;
			self.deviatoricerror_groupthreshold = 0;
			self.deviatoricerror_maximum			= 0;
			self.restart								= 0;
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			md = checkfield(md,'fieldname','amr.level_max','numel',[1],'>=',0,'<=',5);
   		md = checkfield(md,'fieldname','amr.gradation','numel',[1],'>=',1.1,'<=',5.0,'NaN',1);
   		md = checkfield(md,'fieldname','amr.lag','numel',[1],'>=',1.0,'<=',3.0,'NaN',1);
         md = checkfield(md,'fieldname','amr.groundingline_distance','numel',[1],'>=',0,'NaN',1,'Inf',1);
         md = checkfield(md,'fieldname','amr.icefront_distance','numel',[1],'>=',0,'NaN',1,'Inf',1);
         md = checkfield(md,'fieldname','amr.thicknesserror_threshold','numel',[1],'>=',0,'<=',1,'NaN',1);
         md = checkfield(md,'fieldname','amr.thicknesserror_groupthreshold','numel',[1],'>=',0,'<=',1,'NaN',1);
         md = checkfield(md,'fieldname','amr.thicknesserror_maximum','numel',[1],'>=',0,'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','amr.deviatoricerror_threshold','numel',[1],'>=',0,'<=',1,'NaN',1);			
			md = checkfield(md,'fieldname','amr.deviatoricerror_groupthreshold','numel',[1],'>=',0,'<=',1,'NaN',1);			
         md = checkfield(md,'fieldname','amr.deviatoricerror_maximum','numel',[1],'>=',0,'NaN',1,'Inf',1);		
		   md = checkfield(md,'fieldname','amr.restart','numel',[1],'>=',0,'<=',1,'NaN',1);
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   amrneopz parameters:'));

			fielddisplay(self,'level_max',['maximum refinement level (1, 2, 3, 4 or 5)']);
			fielddisplay(self,'gradation',['maximum ratio between two adjacent edges']);
			fielddisplay(self,'lag',['lag used to unrefine the elements']);
         fielddisplay(self,'groundingline_distance',['distance around the grounding line which elements will be refined']);
         fielddisplay(self,'icefront_distance',['distance around the ice front which elements will be refined']);
         fielddisplay(self,'thicknesserror_threshold',['maximum threshold thickness error permitted']);
         fielddisplay(self,'thicknesserror_groupthreshold',['maximum group threshold thickness error permitted']);
         fielddisplay(self,'thicknesserror_maximum',['maximum thickness error permitted']);
			fielddisplay(self,'deviatoricerror_threshold',['maximum threshold deviatoricstress error permitted']);
			fielddisplay(self,'deviatoricerror_groupthreshold',['maximum group threshold deviatoricstress error permitted']);
			fielddisplay(self,'deviatoricerror_maximum',['maximum deviatoricstress error permitted']);
         fielddisplay(self,'restart',['indicates if ReMesh() will call before first time step']);
		end % }}}
		function marshall(self,prefix,md,fid) % {{{

			WriteData(fid,prefix,'name','md.amr.type','data',2,'format','Integer');
			WriteData(fid,prefix,'object',self,'class','amr','fieldname','level_max','format','Integer');
			WriteData(fid,prefix,'object',self,'class','amr','fieldname','gradation','format','Double');
			WriteData(fid,prefix,'object',self,'class','amr','fieldname','lag','format','Double');
         WriteData(fid,prefix,'object',self,'class','amr','fieldname','groundingline_distance','format','Double');
         WriteData(fid,prefix,'object',self,'class','amr','fieldname','icefront_distance','format','Double');
         WriteData(fid,prefix,'object',self,'class','amr','fieldname','thicknesserror_threshold','format','Double');
         WriteData(fid,prefix,'object',self,'class','amr','fieldname','thicknesserror_groupthreshold','format','Double');
         WriteData(fid,prefix,'object',self,'class','amr','fieldname','thicknesserror_maximum','format','Double');
			WriteData(fid,prefix,'object',self,'class','amr','fieldname','deviatoricerror_threshold','format','Double');
			WriteData(fid,prefix,'object',self,'class','amr','fieldname','deviatoricerror_groupthreshold','format','Double');
			WriteData(fid,prefix,'object',self,'class','amr','fieldname','deviatoricerror_maximum','format','Double');
		   WriteData(fid,prefix,'object',self,'class','amr','fieldname','restart','format','Integer');
		end % }}}
		function savemodeljs(self,fid,modelname) % {{{
		
			writejsdouble(fid,[modelname '.amr.level_max'],self.level_max);
			writejsdouble(fid,[modelname '.amr.gradation'],self.gradation);
			writejsdouble(fid,[modelname '.amr.lag'],self.lag);
			writejsdouble(fid,[modelname '.amr.groundingline_distance'],self.groundingline_distance);
			writejsdouble(fid,[modelname '.amr.icefront_distance'],self.icefront_distance);
			writejsdouble(fid,[modelname '.amr.thicknesserror_threshold'],self.thicknesserror_threshold);
			writejsdouble(fid,[modelname '.amr.thicknesserror_groupthreshold'],self.thicknesserror_threshold);
			writejsdouble(fid,[modelname '.amr.thicknesserror_maximum'],self.thicknesserror_maximum);
			writejsdouble(fid,[modelname '.amr.deviatoricerror_threshold'],self.deviatoricerror_threshold);
			writejsdouble(fid,[modelname '.amr.deviatoricerror_groupthreshold'],self.deviatoricerror_threshold);
			writejsdouble(fid,[modelname '.amr.deviatoricerror_maximum'],self.deviatoricerror_maximum);
			writejsdouble(fid,[modelname '.amr.restart'],self.restart);

		end % }}}
	end
end
