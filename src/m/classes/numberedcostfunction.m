%NUMBEREDCOSTFUNCTION class definition
%
%   Usage:
%      numberedcostfunction=numberedcostfunction();

classdef numberedcostfunction
	properties (SetAccess=public) 
		name								 = '';
		definitionstring				 = '';
		cost_functions              = NaN
		cost_functions_coefficients = NaN
		vx_obs                      = NaN
		vy_obs                      = NaN
		vz_obs                      = NaN
		vel_obs                     = NaN
		thickness_obs               = NaN
		surface_obs                 = NaN
	end
	methods
		function self = numberedcostfunction(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				case 1
					self=structtoobj(numberedcostfunction(),varargin{1});
				otherwise
					error('construtor not supported yet');
			end
		end % }}}
		function self = extrude(self,md) % {{{
			self.vx_obs=project3d(md,'vector',self.vx_obs,'type','node');
			self.vy_obs=project3d(md,'vector',self.vy_obs,'type','node');
			self.vel_obs=project3d(md,'vector',self.vel_obs,'type','node');
			self.thickness_obs=project3d(md,'vector',self.thickness_obs,'type','node');
			if numel(self.cost_functions_coefficients)>1,self.cost_functions_coefficients=project3d(md,'vector',self.cost_functions_coefficients,'type','node');end;
		end % }}}
		function self = setdefaultparameters(self) % {{{

			%several responses can be used:
			self.cost_functions=101;
			

		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			num_costfunc=size(self.cost_functions,2);

			md = checkfield(md,'fieldname','md.outputdefinition.definition{X}.dependent{X}.cost_functions','field',self.cost_functions,'size',[1 num_costfunc],'values',supportedcostfunctions());
			md = checkfield(md,'fieldname','md.outputdefinition.definition{X}.cost_functions_coefficients','field',self.cost_functions_coefficients,'size',[md.mesh.numberofvertices numel(self.cost_functions)],'>=',0);
			
			if strcmp(solution,'BalancethicknessSolution')
				md = checkfield(md,'fieldname','md.outputdefinition.definition{X}.thickness_obs','field',self.thickness_obs,'size',[md.mesh.numberofvertices 1],'NaN',1,'Inf',1);
				md = checkfield(md,'fieldname','md.outputdefinition.definition{X}.surface_obs','field',self.surface_obs,'size',[md.mesh.numberofvertices 1],'NaN',1,'Inf',1);
			elseif strcmp(solution,'BalancethicknessSoftSolution')
				md = checkfield(md,'fieldname','md.outputdefinition.definition{X}.thickness_obs','field',self.thickness_obs,'size',[md.mesh.numberofvertices 1],'NaN',1,'Inf',1);
			else
				md = checkfield(md,'fieldname','md.outputdefinition.definition{X}.vx_obs','field',self.vx_obs,'size',[md.mesh.numberofvertices 1],'NaN',1,'Inf',1);
				if ~strcmp(domaintype(md.mesh),'2Dvertical'),
					md = checkfield(md,'fieldname','md.outputdefinition.definition{X}.vy_obs','field',self.vy_obs,'size',[md.mesh.numberofvertices 1],'NaN',1,'Inf',1);
				end
			end
		end % }}}
		function disp(self) % {{{
			disp(sprintf('   numberedcostfunction parameters:'));
			fielddisplay(self,'cost_functions','indicate the type of response for each optimization step');
			fielddisplay(self,'cost_functions_coefficients','cost_functions_coefficients applied to the misfit of each vertex and for each control_parameter');
			fielddisplay(self,'vx_obs','observed velocity x component [m/yr]');
			fielddisplay(self,'vy_obs','observed velocity y component [m/yr]');
			fielddisplay(self,'vel_obs','observed velocity magnitude [m/yr]');
			fielddisplay(self,'thickness_obs','observed thickness [m]');
			fielddisplay(self,'surface_obs','observed surface elevation [m]');
			disp('Available cost functions:');
			disp('   101: SurfaceAbsVelMisfit');
			disp('   102: SurfaceRelVelMisfit');
			disp('   103: SurfaceLogVelMisfit');
			disp('   104: SurfaceLogVxVyMisfit');
			disp('   105: SurfaceAverageVelMisfit');
			disp('   201: ThicknessAbsMisfit');
			disp('   501: DragCoefficientAbsGradient');
			disp('   502: RheologyBbarAbsGradient');
			disp('   503: ThicknessAbsGradient');
		end % }}}
		function marshall(self,prefix,md,fid) % {{{

			yts=md.constants.yts;

			WriteData(fid,prefix,'data',self.name,'name','md.numberedcostfunction.name','format','String');
			WriteData(fid,prefix,'data',self.definitionstring,'name','md.numberedcostfunction.definitionstring','format','String');

			WriteData(fid,prefix,'data',self.cost_functions_coefficients,'name','md.numberedcostfunction.cost_functions_coefficients','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'data',self.vx_obs,'name','md.numberedcostfunction.vx_obs','format','DoubleMat','mattype',1,'scale',1./yts);
			WriteData(fid,prefix,'data',self.vy_obs,'name','md.numberedcostfunction.vy_obs','format','DoubleMat','mattype',1,'scale',1./yts);
			WriteData(fid,prefix,'data',self.vz_obs,'name','md.numberedcostfunction.vz_obs','format','DoubleMat','mattype',1,'scale',1./yts);
			WriteData(fid,prefix,'data',self.vel_obs,'name','md.numberedcostfunction.vel_obs','format','DoubleMat','mattype',1,'scale',1./yts);
			if(numel(self.thickness_obs)==md.mesh.numberofelements),
				mattype=2;
			else
				mattype=1;
			end
			WriteData(fid,prefix,'data',self.thickness_obs,'name','md.numberedcostfunction.thickness_obs','format','DoubleMat','mattype',mattype);
			WriteData(fid,prefix,'data',self.surface_obs,'name','md.numberedcostfunction.surface_obs','format','DoubleMat','mattype',mattype);

			%process cost functions
			num_cost_functions=size(self.cost_functions,2);
			data=marshallcostfunctions(self.cost_functions);
			WriteData(fid,prefix,'data',data,'name','md.numberedcostfunction.cost_functions','format','StringArray');
			WriteData(fid,prefix,'data',num_cost_functions,'name','md.numberedcostfunction.num_cost_functions','format','Integer');
		end % }}}
		function savemodeljs(self,fid,modelname) % {{{
		
			writejs2Darray(fid,[modelname '.inversion.cost_functions'],self.cost_functions);
			writejs2Darray(fid,[modelname '.inversion.cost_functions_coefficients'],self.cost_functions_coefficients);
			writejs1Darray(fid,[modelname '.inversion.vx_obs'],self.vx_obs);
			writejs1Darray(fid,[modelname '.inversion.vy_obs'],self.vy_obs);
			writejs1Darray(fid,[modelname '.inversion.vz_obs'],self.vz_obs);
			writejs1Darray(fid,[modelname '.inversion.vel_obs'],self.vel_obs);
			writejs1Darray(fid,[modelname '.inversion.thickness_obs'],self.thickness_obs);
			writejs1Darray(fid,[modelname '.inversion.surface_obs'],self.surface_obs);

		end % }}}
	end
end
