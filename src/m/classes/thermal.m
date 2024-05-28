%THERMAL class definition
%
%   Usage:
%      thermal=thermal();

classdef thermal
	properties (SetAccess=public) 
		spctemperature    = NaN;
		penalty_threshold = 0;
		stabilization     = 0;
		reltol				= 0;
		maxiter           = 0;
		penalty_lock      = 0;
		penalty_factor    = 0;
		isenthalpy        = 0;
		isdynamicbasalspc = 0;
		isdrainicecolumn   = 0;
		watercolumn_upperlimit= 0;
		fe                = 'P1';
		requested_outputs = {};
	end
	methods
		function self = extrude(self,md) % {{{
			self.spctemperature=project3d(md,'vector',self.spctemperature,'type','node','layer',md.mesh.numberoflayers,'padding',NaN);
			if (length(md.initialization.temperature)==md.mesh.numberofvertices),
				self.spctemperature=NaN(md.mesh.numberofvertices,1);
				pos=find(md.mesh.vertexonsurface);
				self.spctemperature(pos)=md.initialization.temperature(pos); %impose observed temperature on surface
			end
		end % }}}
		function self = thermal(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function list = defaultoutputs(self,md) % {{{

			if self.isenthalpy,
				list = {'Enthalpy','Temperature','Waterfraction','Watercolumn','BasalforcingsGroundediceMeltingRate'};
			else
				list = {'Temperature','BasalforcingsGroundediceMeltingRate'};
			end

		end % }}}
		function self = setdefaultparameters(self) % {{{

			%Number of unstable constraints acceptable
			self.penalty_threshold=0;

			%Type of stabilization used
			self.stabilization=1;

			%Relative tolerance for the enthalpy convergence
			self.reltol=0.01;

			%Maximum number of iterations
			self.maxiter=100;

			%factor used to compute the values of the penalties: kappa=max(stiffness matrix)*10^penalty_factor
			self.penalty_factor=3;

			%Should we use cold ice (default) or enthalpy formulation
			self.isenthalpy=0;

			%will basal boundary conditions be set dynamically
			self.isdynamicbasalspc=0;
		
			%wether waterfraction drainage is enabled
			self.isdrainicecolumn=1;

			%set an upper limit for local stored watercolumn
			self.watercolumn_upperlimit=1000;

			%Linear elements by default
			self.fe='P1';

			%default output
			self.requested_outputs={'default'};
		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			%Early return
			if (~ismember('ThermalAnalysis',analyses) & ~ismember('EnthalpyAnalysis',analyses)) | (strcmp(solution,'TransientSolution') & md.transient.isthermal==0), return; end

			md = checkfield(md,'fieldname','thermal.stabilization','numel',[1],'values',[0 1 2 3]);
			md = checkfield(md,'fieldname','thermal.spctemperature','Inf',1,'timeseries',1,'>=',0);
			md = checkfield(md,'fieldname','thermal.fe','values',{'P1','P1xP2','P1xP3'});
			if (ismember('EnthalpyAnalysis',analyses) & md.thermal.isenthalpy & dimension(md.mesh)==3),
				md = checkfield(md,'fieldname','thermal.isdrainicecolumn','numel',[1],'values',[0 1]);
				md = checkfield(md,'fieldname','thermal.watercolumn_upperlimit','>=',0);

				%Make sure the spc are less than melting point (Josh commented out the next 2 lines)
				TEMP=md.thermal.spctemperature(1:md.mesh.numberofvertices,:);
				replicate=repmat(md.geometry.surface-md.mesh.z,1,size(md.thermal.spctemperature,2));
				pos=find(~isnan(TEMP));
				md = checkfield(md,'fieldname','thermal.spctemperature','field',TEMP(pos)',...
					'<=',md.materials.meltingpoint-md.materials.beta*md.materials.rho_ice*md.constants.g*replicate(pos)+10^-5,...
					'message','spctemperature should be less or equal than the adjusted melting point');

				md = checkfield(md,'fieldname','thermal.isenthalpy','numel',[1],'values',[0 1]);
				md = checkfield(md,'fieldname','thermal.isdynamicbasalspc','numel', [1],'values',[0 1]);
				if(md.thermal.isenthalpy)
					if isnan(md.stressbalance.reltol),
						md = checkmessage(md,['for a steadystate computation, thermal.reltol (relative convergence criterion) must be defined!']);
					end 
					md = checkfield(md,'fieldname','thermal.reltol','>',0.,'message','reltol must be larger than zero');
				end
			end

		 md = checkfield(md,'fieldname','thermal.requested_outputs','stringrow',1);
    end % }}} 
		function disp(self) % {{{
			disp(sprintf('   Thermal solution parameters:'));

			fielddisplay(self,'spctemperature','temperature constraints (NaN means no constraint) [K]');
			fielddisplay(self,'stabilization','0: no, 1: artificial_diffusivity, 2: SUPG, 3: anisotropic SUPG');
			fielddisplay(self,'reltol','relative tolerance convergence criterion for enthalpy');
			fielddisplay(self,'maxiter','maximum number of non linear iterations');
			fielddisplay(self,'penalty_lock','stabilize unstable thermal constraints that keep zigzagging after n iteration (default is 0, no stabilization)');
			fielddisplay(self,'penalty_threshold','threshold to declare convergence of thermal solution (default is 0)');
			fielddisplay(self,'penalty_factor','scaling exponent (default is 3)');
			fielddisplay(self,'isenthalpy','use an enthalpy formulation to include temperate ice (default is 0)');
			fielddisplay(self,'isdynamicbasalspc','enable dynamic setting of basal forcing. required for enthalpy formulation (default is 0)');
			fielddisplay(self,'isdrainicecolumn','wether waterfraction drainage is enabled for enthalpy formulation (default is 1)'); 
			fielddisplay(self,'watercolumn_upperlimit','upper limit of basal watercolumn for enthalpy formulation (default is 1000m)');
			fielddisplay(self,'fe','Finite Element type: ''P1'' (default), ''P1xP2''');
			fielddisplay(self,'requested_outputs','additional outputs requested');

		end % }}}
		function marshall(self,prefix,md,fid) % {{{
			WriteData(fid,prefix,'object',self,'fieldname','spctemperature','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'fieldname','penalty_threshold','format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','stabilization','format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','reltol','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','maxiter','format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','penalty_lock','format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','penalty_factor','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','isenthalpy','format','Boolean');
			WriteData(fid,prefix,'object',self,'fieldname','isdrainicecolumn','format','Boolean');
			WriteData(fid,prefix,'object',self,'fieldname','watercolumn_upperlimit','format','Double');
			WriteData(fid,prefix,'object',self,'fieldname','fe','format','String');
			WriteData(fid,prefix,'object',self,'fieldname','isdynamicbasalspc','format','Boolean');

			%process requested outputs
			outputs = self.requested_outputs;
			pos  = find(ismember(outputs,'default'));
			if ~isempty(pos),
				outputs(pos) = [];                         %remove 'default' from outputs
				outputs      = [outputs defaultoutputs(self,md)]; %add defaults
			end
			WriteData(fid,prefix,'data',outputs,'name','md.thermal.requested_outputs','format','StringArray');
        	end % }}}
		function savemodeljs(self,fid,modelname) % {{{
		
			writejs1Darray(fid,[modelname '.thermal.spctemperature'],self.spctemperature);
			writejsdouble(fid,[modelname '.thermal.penalty_threshold'],self.penalty_threshold);
			writejsdouble(fid,[modelname '.thermal.stabilization'],self.stabilization);
			writejsdouble(fid,[modelname '.thermal.reltol'],self.reltol);
			writejsdouble(fid,[modelname '.thermal.maxiter'],self.maxiter);
			writejsdouble(fid,[modelname '.thermal.penalty_lock'],self.penalty_lock);
			writejsdouble(fid,[modelname '.thermal.penalty_factor'],self.penalty_factor);
			writejsdouble(fid,[modelname '.thermal.isenthalpy'],self.isenthalpy);
			writejsdouble(fid,[modelname '.thermal.isdrainicecolumn'],self.isdrainicecolumn);
			writejsdouble(fid,[modelname '.thermal.watercolumn_upperlimit'],self.watercolumn_upperlimit);
			writejsdouble(fid,[modelname '.thermal.isdynamicbasalspc'],self.isdynamicbasalspc);
			writejscellstring(fid,[modelname '.thermal.requested_outputs'],self.requested_outputs);

		end % }}}
	end
end
