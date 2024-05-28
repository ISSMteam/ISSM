%AMR Class definition
%
%   Usage:
%      md.amr=amr();

classdef amr
	properties (SetAccess=public) 
		hmin = 0.; 
		hmax = 0.;
		fieldname = '';
		err = 0.;
		keepmetric = 0;
		gradation = 0.;
		groundingline_resolution = 0.;
		groundingline_distance = 0.;
		icefront_resolution = 0.;
		icefront_distance = 0.;
		thicknesserror_resolution = 0.;
		thicknesserror_threshold = 0.;
		thicknesserror_groupthreshold = 0.;
		thicknesserror_maximum = 0.;
		deviatoricerror_resolution = 0.;
		deviatoricerror_threshold = 0.;
		deviatoricerror_groupthreshold = 0.;
		deviatoricerror_maximum = 0.;
		restart=0.;
	end
	methods (Static)
 		function self = loadobj(self) % {{{
         % This function is directly called by matlab when a model object is
         % loaded. Update old properties here

         if verLessThan('matlab','7.9'),
            disp('Warning: your matlab version is old and there is a risk that load does not work correctly');
            disp('         if the model is not loaded correctly, rename temporarily loadobj so that matlab does not use it');

            % This is a Matlab bug: all the fields of md have their default value
            % Example of error message:
            % Warning: Error loading an object of class 'model':
            % Undefined function or method 'exist' for input arguments of type 'cell'
            %
            % This has been fixed in MATLAB 7.9 (R2009b) and later versions
         end

         %2017 September 15th
         if isstruct(self),
            disp('WARNING: updating amr. Now the default is amr with bamg');
            disp('         some old fields were not converted');
            disp('         see the new fields typing md.amr and modify them properly');
            obj2 = self;
            self = amr();
            %Converting region_level_max to groundingline_distance
            if(obj2.region_level_max>0 && obj2.region_level_1>obj2.region_level_max)
               self.groundingline_distance	= obj2.region_level_max;
					self.keepmetric					= 0;
					self.fieldname						= 'None';
            end
         end
      end% }}}
	end
	methods
		function self = amr(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = setdefaultparameters(self) % {{{

			%hmin and hmax
			self.hmin=100.;
			self.hmax=100.e3;

			%fields
			self.fieldname ='Vel';
			self.err=3.;

			%keep metric?
			self.keepmetric=1;

			%control of element lengths
			self.gradation=1.5;

			%other criteria
			self.groundingline_resolution=500.;
			self.groundingline_distance=0.;
			self.icefront_resolution=500.;
			self.icefront_distance=0.;
			self.thicknesserror_resolution=500.;
			self.thicknesserror_threshold=0.;
			self.thicknesserror_groupthreshold=0.;
			self.thicknesserror_maximum=0.;
			self.deviatoricerror_resolution=500.;
			self.deviatoricerror_threshold=0.;
			self.deviatoricerror_groupthreshold=0.;
			self.deviatoricerror_maximum=0.;
			
			%is restart? This calls femmodel->ReMesh() before first time step. 
			self.restart=0;

		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			md = checkfield(md,'fieldname','amr.hmax','numel',[1],'>',0,'NaN',1);
			md = checkfield(md,'fieldname','amr.hmin','numel',[1],'>',0,'<',self.hmax,'NaN',1);
			%md = checkfield(md,'fieldname','amr.fieldname','string',[1]);
			md = checkfield(md,'fieldname','amr.keepmetric','numel',[1],'>=',0,'<=',1,'NaN',1);
			md = checkfield(md,'fieldname','amr.gradation','numel',[1],'>=',1.1,'<=',5,'NaN',1);
			md = checkfield(md,'fieldname','amr.groundingline_resolution','numel',[1],'>',0,'<',self.hmax,'NaN',1);
			md = checkfield(md,'fieldname','amr.groundingline_distance','numel',[1],'>=',0,'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','amr.icefront_resolution','numel',[1],'>',0,'<',self.hmax,'NaN',1);
			md = checkfield(md,'fieldname','amr.icefront_distance','numel',[1],'>=',0,'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','amr.thicknesserror_resolution','numel',[1],'>',0,'<',self.hmax,'NaN',1);
			md = checkfield(md,'fieldname','amr.thicknesserror_threshold','numel',[1],'>=',0,'<=',1,'NaN',1);
			md = checkfield(md,'fieldname','amr.thicknesserror_groupthreshold','numel',[1],'>=',0,'<=',1,'NaN',1);
			md = checkfield(md,'fieldname','amr.thicknesserror_maximum','numel',[1],'>=',0,'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','amr.deviatoricerror_resolution','numel',[1],'>',0,'<',self.hmax,'NaN',1);
			md = checkfield(md,'fieldname','amr.deviatoricerror_threshold','numel',[1],'>=',0,'<=',1,'NaN',1);
			md = checkfield(md,'fieldname','amr.deviatoricerror_groupthreshold','numel',[1],'>=',0,'<=',1,'NaN',1);
			md = checkfield(md,'fieldname','amr.deviatoricerror_maximum','numel',[1],'>=',0,'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','amr.restart','numel',[1],'>=',0,'<=',1,'NaN',1);
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   amr parameters:'));

			fielddisplay(self,'hmin',['minimum element length']);
			fielddisplay(self,'hmax',['maximum element length']);
			fielddisplay(self,'fieldname',['name of input that will be used to compute the metric (should be an input of FemModel)']);
			fielddisplay(self,'keepmetric',['indicates whether the metric should be kept every remeshing time']);
			fielddisplay(self,'gradation',['maximum ratio between two adjacent edges']);
			fielddisplay(self,'groundingline_resolution',['element length near the grounding line']);
			fielddisplay(self,'groundingline_distance',['distance around the grounding line which elements will be refined']);
			fielddisplay(self,'icefront_resolution',['element length near the ice front']);
			fielddisplay(self,'icefront_distance',['distance around the ice front which elements will be refined']);
			fielddisplay(self,'thicknesserror_resolution',['element length when thickness error estimator is used']);
			fielddisplay(self,'thicknesserror_threshold',['maximum threshold thickness error permitted']);
			fielddisplay(self,'thicknesserror_groupthreshold',['maximum group threshold thickness error permitted']);
			fielddisplay(self,'thicknesserror_maximum',['maximum thickness error permitted']);
			fielddisplay(self,'deviatoricerror_resolution',['element length when deviatoric stress error estimator is used']);
			fielddisplay(self,'deviatoricerror_threshold',['maximum threshold deviatoricstress error permitted']);
			fielddisplay(self,'deviatoricerror_groupthreshold',['maximum group threshold deviatoricstress error permitted']);
			fielddisplay(self,'deviatoricerror_maximum',['maximum deviatoricstress error permitted']);
			fielddisplay(self,'restart',['indicates if ReMesh() will call before first time step']);
		end % }}}
			function savemodeljs(self,fid,modelname) % {{{
			%error('not implemented yet!');
		end % }}}
	function marshall(self,prefix,md,fid) % {{{

			WriteData(fid,prefix,'name','md.amr.type','data',1,'format','Integer');
			WriteData(fid,prefix,'object',self,'class','amr','fieldname','hmin','format','Double');
			WriteData(fid,prefix,'object',self,'class','amr','fieldname','hmax','format','Double');
			WriteData(fid,prefix,'object',self,'class','amr','fieldname','fieldname','format','String');
			WriteData(fid,prefix,'object',self,'class','amr','fieldname','err','format','Double');
			WriteData(fid,prefix,'object',self,'class','amr','fieldname','keepmetric','format','Integer');
			WriteData(fid,prefix,'object',self,'class','amr','fieldname','gradation','format','Double');
			WriteData(fid,prefix,'object',self,'class','amr','fieldname','groundingline_resolution','format','Double');
			WriteData(fid,prefix,'object',self,'class','amr','fieldname','groundingline_distance','format','Double');
			WriteData(fid,prefix,'object',self,'class','amr','fieldname','icefront_resolution','format','Double');
			WriteData(fid,prefix,'object',self,'class','amr','fieldname','icefront_distance','format','Double');
			WriteData(fid,prefix,'object',self,'class','amr','fieldname','thicknesserror_resolution','format','Double');
			WriteData(fid,prefix,'object',self,'class','amr','fieldname','thicknesserror_threshold','format','Double');
			WriteData(fid,prefix,'object',self,'class','amr','fieldname','thicknesserror_groupthreshold','format','Double');
			WriteData(fid,prefix,'object',self,'class','amr','fieldname','thicknesserror_maximum','format','Double');
			WriteData(fid,prefix,'object',self,'class','amr','fieldname','deviatoricerror_resolution','format','Double');
			WriteData(fid,prefix,'object',self,'class','amr','fieldname','deviatoricerror_threshold','format','Double');
			WriteData(fid,prefix,'object',self,'class','amr','fieldname','deviatoricerror_groupthreshold','format','Double');
			WriteData(fid,prefix,'object',self,'class','amr','fieldname','deviatoricerror_maximum','format','Double');
			WriteData(fid,prefix,'object',self,'class','amr','fieldname','restart','format','Integer');

		end % }}}
	end
end
