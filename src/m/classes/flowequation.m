%FLOWEQUATION class definition
%
%   Usage:
%      flowequation=flowequation();

classdef flowequation
	properties (SetAccess=public) 
		isSIA                          = 0;
		isSSA                          = 0;
		isL1L2                         = 0;
		isMOLHO                         = 0;
		isHO                           = 0;
		isFS                           = 0;
		isNitscheBC                    = 0;
		FSNitscheGamma                 = 1e6;
		fe_SSA                         = '';
		fe_HO                          = '';
		fe_FS                          = '';
		augmented_lagrangian_r         = 1.;
		augmented_lagrangian_rhop      = 1.;
		augmented_lagrangian_rlambda   = 1.;
		augmented_lagrangian_rholambda = 1.;
		XTH_theta                      = 0.;
		vertex_equation                = NaN;
		element_equation               = NaN;
		borderSSA                      = NaN;
		borderHO                       = NaN;
		borderFS                       = NaN;
	end
	methods (Static)
		function self = loadobj(self) % {{{
			% This function is directly called by matlab when a model object is
			% loaded. If the input is a struct it is an old version of this class and
			% old fields must be recovered (make sure they are in the deprecated
			% model properties)

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

			if isstruct(self)
				disp('Recovering flowequation from older version');
				objstruct = self;
				self = structtoobj(flowequation(),objstruct);

				%2013 July 23rd
				if isfield(objstruct,'ishutter'),      self.isSIA     = objstruct.ishutter;       end; 
				if isfield(objstruct,'ismacayeal'),    self.isSSA     = objstruct.ismacayeal;     end; 
				if isfield(objstruct,'ispattyn'),      self.isHO      = objstruct.ispattyn;       end; 
				if isfield(objstruct,'isstokes'),      self.isFS      = objstruct.isstokes;       end; 
				if isfield(objstruct,'bordermacayeal'),self.borderSSA = objstruct.bordermacayeal; end; 
				if isfield(objstruct,'borderpattyn'),  self.borderHO  = objstruct.borderpattyn;   end; 
				if isfield(objstruct,'borderstokes'),  self.borderFS  = objstruct.borderstokes;   end; 

				%May 31 2022
				if isfield(objstruct,'isMLHO')
					self.isMOLHO = objstruct.isMLHO;
				end
			end


		end% }}}
	end
	methods
		function self = extrude(self,md) % {{{
			self.element_equation=project3d(md,'vector',self.element_equation,'type','element');
			self.vertex_equation=project3d(md,'vector',self.vertex_equation,'type','node');
			self.borderSSA=project3d(md,'vector',self.borderSSA,'type','node');
			self.borderHO=project3d(md,'vector',self.borderHO,'type','node');
			self.borderFS=project3d(md,'vector',self.borderFS,'type','node');
		end % }}}
		function self = flowequation(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = setdefaultparameters(self) % {{{

			%P1 for SSA
			self.fe_SSA= 'P1';

			%P1 for HO
			self.fe_HO= 'P1';

			%MINI condensed element for FS by default
			self.fe_FS = 'MINIcondensed';
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			%Early return
			if (~ismember('StressbalanceAnalysis',analyses) & ~ismember('StressbalanceSIAAnalysis',analyses)) | (strcmp(solution,'TransientSolution') & md.transient.isstressbalance==0), return; end

			md = checkfield(md,'fieldname','flowequation.isSIA','numel',[1],'values',[0 1]);
			md = checkfield(md,'fieldname','flowequation.isSSA','numel',[1],'values',[0 1]);
			md = checkfield(md,'fieldname','flowequation.isL1L2','numel',[1],'values',[0 1]);
			md = checkfield(md,'fieldname','flowequation.isMOLHO','numel',[1],'values',[0 1]);
			md = checkfield(md,'fieldname','flowequation.isHO','numel',[1],'values',[0 1]);
			md = checkfield(md,'fieldname','flowequation.isFS','numel',[1],'values',[0 1]);
			md = checkfield(md,'fieldname','flowequation.isNitscheBC','numel',[1],'values',[0 1]);
			md = checkfield(md,'fieldname','flowequation.FSNitscheGamma','numel',[1], '>=', 0.);
			md = checkfield(md,'fieldname','flowequation.fe_SSA','values',{'P1','P1bubble','P1bubblecondensed','P2','P2bubble'});
			md = checkfield(md,'fieldname','flowequation.fe_HO' ,'values',{'P1','P1bubble','P1bubblecondensed','P1xP2','P2xP1','P2','P2bubble','P1xP3','P1xP4','P2xP4'});
			md = checkfield(md,'fieldname','flowequation.fe_FS' ,'values',{'P1P1','P1P1GLS','MINIcondensed','MINI','TaylorHood','LATaylorHood','XTaylorHood','OneLayerP4z','CrouzeixRaviart','LACrouzeixRaviart'});
			md = checkfield(md,'fieldname','flowequation.augmented_lagrangian_r','numel',[1],'>=',0.);
			md = checkfield(md,'fieldname','flowequation.augmented_lagrangian_rlambda','numel',[1],'>=',0.);
			md = checkfield(md,'fieldname','flowequation.augmented_lagrangian_rhop','numel',[1],'>=',0.);
			md = checkfield(md,'fieldname','flowequation.augmented_lagrangian_rholambda','numel',[1],'>=',0.);
			md = checkfield(md,'fieldname','flowequation.XTH_theta','numel',[1],'>=',0.,'<',0.5);
			md = checkfield(md,'fieldname','flowequation.borderSSA','size',[md.mesh.numberofvertices 1],'values',[0 1]);
			md = checkfield(md,'fieldname','flowequation.borderHO','size',[md.mesh.numberofvertices 1],'values',[0 1]);
			md = checkfield(md,'fieldname','flowequation.borderFS','size',[md.mesh.numberofvertices 1],'values',[0 1]);
			if strcmp(domaintype(md.mesh),'2Dhorizontal')
				md = checkfield(md,'fieldname','flowequation.vertex_equation','size',[md.mesh.numberofvertices 1],'values',[1,2,4]);
				md = checkfield(md,'fieldname','flowequation.element_equation','size',[md.mesh.numberofelements 1],'values',[1,2,4]);
			elseif strcmp(domaintype(md.mesh),'3Dsurface')
				md = checkfield(md,'fieldname','flowequation.vertex_equation','size',[md.mesh.numberofvertices 1],'values',[1:2]);
				md = checkfield(md,'fieldname','flowequation.element_equation','size',[md.mesh.numberofelements 1],'values',[1:2]);
			elseif strcmp(domaintype(md.mesh),'2Dvertical')
				md = checkfield(md,'fieldname','flowequation.vertex_equation','size',[md.mesh.numberofvertices 1],'values',[2,5,6]);
				md = checkfield(md,'fieldname','flowequation.element_equation','size',[md.mesh.numberofelements 1],'values',[2,5,6]);
			elseif strcmp(domaintype(md.mesh),'3D'),
				md = checkfield(md,'fieldname','flowequation.vertex_equation','size',[md.mesh.numberofvertices 1],'values',[0:9]);
				md = checkfield(md,'fieldname','flowequation.element_equation','size',[md.mesh.numberofelements 1],'values',[0:9]);
			else
				error('Case not supported yet');
			end
			if ~(self.isSIA || self.isSSA || self.isL1L2 || self.isMOLHO || self.isHO || self.isFS),
				md = checkmessage(md,['no element types set for this model']);
			end
			if ismember('StressbalanceSIAAnalysis',analyses),
				if any(self.element_equation==1),
					if(self.vertex_equation & md.mask.ocean_levelset<0.),
						disp(sprintf('\n !!! Warning: SIA''s model is not consistent on ice shelves !!!\n'));
					end
				end
			end

		end % }}}
		function disp(self) % {{{
			disp(sprintf('   flow equation parameters:'));

			fielddisplay(self,'isSIA','is the Shallow Ice Approximation (SIA) used?');
			fielddisplay(self,'isSSA','is the Shelfy-Stream Approximation (SSA) used?');
			fielddisplay(self,'isL1L2','is the L1L2 approximation used?');
			fielddisplay(self,'isMOLHO','is the MOno-Layer Higher-Order (MOLHO) approximation used?');
			fielddisplay(self,'isHO','is the Higher-Order (HO) approximation used?');
			fielddisplay(self,'isFS','are the Full-FS (FS) equations used?');
			fielddisplay(self,'isNitscheBC','is weakly imposed condition used?');
			fielddisplay(self,'FSNitscheGamma','Gamma value for the Nitsche term (default: 1e6)');
			fielddisplay(self,'fe_SSA','Finite Element for SSA  ''P1'', ''P1bubble'' ''P1bubblecondensed'' ''P2''');
			fielddisplay(self,'fe_HO', 'Finite Element for HO   ''P1'' ''P1bubble'' ''P1bubblecondensed'' ''P1xP2'' ''P2xP1'' ''P2''');
			fielddisplay(self,'fe_FS', 'Finite Element for FS   ''P1P1'' (debugging only) ''P1P1GLS'' ''MINIcondensed'' ''MINI'' ''TaylorHood'' ''XTaylorHood''');
			fielddisplay(self,'vertex_equation','flow equation for each vertex');
			fielddisplay(self,'element_equation','flow equation for each element');
			fielddisplay(self,'borderSSA','vertices on SSA''s border (for tiling)');
			fielddisplay(self,'borderHO','vertices on HO''s border (for tiling)');
			fielddisplay(self,'borderFS','vertices on FS'' border (for tiling)');

		end % }}}
		function marshall(self,prefix,md,fid) % {{{
			WriteData(fid,prefix,'object',self,'fieldname','isSIA','format','Boolean');
			WriteData(fid,prefix,'object',self,'fieldname','isSSA','format','Boolean');
			WriteData(fid,prefix,'object',self,'fieldname','isL1L2','format','Boolean');
			WriteData(fid,prefix,'object',self,'fieldname','isMOLHO','format','Boolean');
			WriteData(fid,prefix,'object',self,'fieldname','isHO','format','Boolean');
			WriteData(fid,prefix,'object',self,'fieldname','isFS','format','Boolean');
			WriteData(fid,prefix,'object',self,'fieldname','isNitscheBC','format','Boolean');
			WriteData(fid,prefix,'object',self,'fieldname','FSNitscheGamma','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','fe_SSA','data',self.fe_SSA,'format','String');
			WriteData(fid,prefix,'object',self,'fieldname','fe_HO' ,'data',self.fe_HO,'format','String');
			WriteData(fid,prefix,'object',self,'fieldname','fe_FS' ,'data',self.fe_FS,'format','String');
			WriteData(fid,prefix,'object',self,'fieldname','augmented_lagrangian_r','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','augmented_lagrangian_rhop','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','augmented_lagrangian_rlambda','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','augmented_lagrangian_rholambda','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','XTH_theta','data',self.XTH_theta ,'format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','borderSSA','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'fieldname','borderHO','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'fieldname','borderFS','format','DoubleMat','mattype',1);
			%convert approximations to enums
			WriteData(fid,prefix,'data',self.vertex_equation,'name','md.flowequation.vertex_equation','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'data',self.element_equation,'name','md.flowequation.element_equation','format','DoubleMat','mattype',2);
		end % }}}
		function savemodeljs(self,fid,modelname) % {{{
		
			writejsdouble(fid,[modelname '.flowequation.isSIA'],self.isSIA);
			writejsdouble(fid,[modelname '.flowequation.isSSA'],self.isSSA);
			writejsdouble(fid,[modelname '.flowequation.isL1L2'],self.isL1L2);
			writejsdouble(fid,[modelname '.flowequation.isMOLHO'],self.isMOLHO);
			writejsdouble(fid,[modelname '.flowequation.isHO'],self.isHO);
			writejsdouble(fid,[modelname '.flowequation.isFS'],self.isFS);
         writejsstring(fid,[modelname '.flowequation.isNitscheBC'],self.isNitscheBC);
         writejsstring(fid,[modelname '.flowequation.FSNitscheGamma'],self.FSNitscheGamma);
         writejsstring(fid,[modelname '.flowequation.fe_SSA'],self.fe_SSA);
			writejsstring(fid,[modelname '.flowequation.fe_HO'],self.fe_HO);
			writejsstring(fid,[modelname '.flowequation.fe_FS'],self.fe_FS);
			writejsdouble(fid,[modelname '.flowequation.augmented_lagrangian_r'],self.augmented_lagrangian_r);
			writejsdouble(fid,[modelname '.flowequation.augmented_lagrangian_rhop'],self.augmented_lagrangian_rhop);
			writejsdouble(fid,[modelname '.flowequation.augmented_lagrangian_rlambda'],self.augmented_lagrangian_rlambda);
			writejsdouble(fid,[modelname '.flowequation.augmented_lagrangian_rholambda'],self.augmented_lagrangian_rholambda);
			writejsdouble(fid,[modelname '.flowequation.XTH_theta'],self.XTH_theta);
			writejs1Darray(fid,[modelname '.flowequation.vertex_equation'],self.vertex_equation);
			writejs1Darray(fid,[modelname '.flowequation.element_equation'],self.element_equation);
			writejs1Darray(fid,[modelname '.flowequation.borderSSA'],self.borderSSA);
			writejs1Darray(fid,[modelname '.flowequation.borderHO'],self.borderHO);
			writejs1Darray(fid,[modelname '.flowequation.borderFS'],self.borderFS);

		end % }}}
	end
end
