%STRESSBALANCE class definition
%
%   Usage:
%      stressbalance=stressbalance();

classdef stressbalance
	properties (SetAccess=public)
		spcvx                  = NaN;
		spcvy                  = NaN;
		spcvx_base             = NaN;
		spcvy_base             = NaN;
		spcvx_shear            = NaN;
		spcvy_shear            = NaN;
		spcvz                  = NaN;
		restol                 = 0;
		reltol                 = 0;
		abstol                 = 0;
		isnewton               = 0;
		FSreconditioning       = 0;
		maxiter                = 0;
		shelf_dampening        = 0;
		vertex_pairing         = NaN;
		penalty_factor         = NaN;
		rift_penalty_lock      = NaN;
		rift_penalty_threshold = 0;
		referential            = NaN;
		loadingforce           = NaN;
		requested_outputs      = {};
	end
	methods
		function self = extrude(self,md) % {{{
			self.spcvx=project3d(md,'vector',self.spcvx,'type','node');
			self.spcvy=project3d(md,'vector',self.spcvy,'type','node');
			self.spcvz=project3d(md,'vector',self.spcvz,'type','node');
			self.referential=project3d(md,'vector',self.referential,'type','node');
			self.loadingforce=project3d(md,'vector',self.loadingforce,'type','node');

			% for MOLHO
			if md.flowequation.isMOLHO
				self.spcvx_base=project3d(md,'vector',self.spcvx_base,'type','node');
				self.spcvy_base=project3d(md,'vector',self.spcvy_base,'type','node');
				self.spcvx_shear=project3d(md,'vector',self.spcvx_shear,'type','poly','degree',4);
				self.spcvy_shear=project3d(md,'vector',self.spcvy_shear,'type','poly','degree',4);
			end
		end % }}}
		function self = stressbalance(varargin) % {{{
			switch nargin
				case 0
					self=setdefaultparameters(self);
				case 1
					inputstruct=varargin{1};
					list1 = properties('stressbalance');
					list2 = fieldnames(inputstruct);
					for i=1:length(list1)
						fieldname = list1{i};
						if ismember(fieldname,list2),
							self.(fieldname) = inputstruct.(fieldname);
						end
					end
				otherwise
					error('constructor not supported');
			end
		end % }}}
		function self = setdefaultparameters(self) % {{{

			%maximum of non-linear iterations.
			self.maxiter=100;

			%Convergence criterion: absolute, relative and residual
			self.restol=10^-4;
			self.reltol=0.01;
			self.abstol=10;

			self.FSreconditioning=10^13;
			self.shelf_dampening=0;

			%Penalty factor applied kappa=max(stiffness matrix)*10^penalty_factor
			self.penalty_factor=3;

			%Stop the iterations of rift if below a threshold
			self.rift_penalty_threshold=0;

			%in some solutions, it might be needed to stop a run when only
			%a few constraints remain unstable. For thermal computation, this
			%parameter is often used.
			self.rift_penalty_lock=10;

			%output default:
			self.requested_outputs={'default'};

		end % }}}
		function md = checkconsistency(self,md,solution,analyses) % {{{

			%Early return
			if ~ismember('StressbalanceAnalysis',analyses), return; end
			if (strcmp(solution,'TransientSolution') & md.transient.isstressbalance == 0), return; end

			md = checkfield(md,'fieldname','stressbalance.spcvx','Inf',1,'timeseries',1);
			md = checkfield(md,'fieldname','stressbalance.spcvy','Inf',1,'timeseries',1);
			md = checkfield(md,'fieldname','stressbalance.spcvz','Inf',1,'timeseries',1);
			md = checkfield(md,'fieldname','stressbalance.restol','size',[1 1],'>',0,'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','stressbalance.reltol','size',[1 1]);
			md = checkfield(md,'fieldname','stressbalance.abstol','size',[1 1]);
			md = checkfield(md,'fieldname','stressbalance.isnewton','numel',[1],'values',[0 1 2]);
			md = checkfield(md,'fieldname','stressbalance.FSreconditioning','size',[1 1],'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','stressbalance.maxiter','size',[1 1],'>=',1);
			md = checkfield(md,'fieldname','stressbalance.referential','size',[md.mesh.numberofvertices 6]);
			md = checkfield(md,'fieldname','stressbalance.loadingforce','size',[md.mesh.numberofvertices 3]);
			md = checkfield(md,'fieldname','stressbalance.requested_outputs','stringrow',1);
			if ~any(isnan(md.stressbalance.vertex_pairing)),
				md = checkfield(md,'fieldname','stressbalance.vertex_pairing','>',0);
			end
			%singular solution
			if ((~(any(~isnan(md.stressbalance.spcvx)) | any(~isnan(md.stressbalance.spcvy)))) & ~any(md.mask.ocean_levelset>0)),
				disp(sprintf('\n !!! Warning: no spc applied, model might not be well posed if no basal friction is applied, check for solution crash\n'));
			end
			%CHECK THAT EACH LINE CONTAINS ONLY NAN VALUES OR NO NAN VALUES
			if any(sum(isnan(md.stressbalance.referential),2)~=0 & sum(isnan(md.stressbalance.referential),2)~=6),
				md = checkmessage(md,['Each line of stressbalance.referential should contain either only NaN values or no NaN values']);
			end
			%CHECK THAT THE TWO VECTORS PROVIDED ARE ORTHOGONAL
			if any(sum(isnan(md.stressbalance.referential),2)==0),
				pos=find(sum(isnan(md.stressbalance.referential),2)==0);
				if any(abs(dot(md.stressbalance.referential(pos,1:3),md.stressbalance.referential(pos,4:6),2))>eps),
					md = checkmessage(md,['Vectors in stressbalance.referential (columns 1 to 3 and 4 to 6) must be orthogonal']);
				end
			end
			%CHECK THAT NO rotation specified for FS Grounded ice at base
			if strcmp(domaintype(md.mesh),'3D') & md.flowequation.isFS,
				pos=find(md.mask.ocean_levelset>0. & md.mesh.vertexonbase);
				if any(~isnan(md.stressbalance.referential(pos,:))),
					md = checkmessage(md,['no referential should be specified for basal vertices of grounded ice']);
				end
				md = checkfield(md,'fieldname','stressbalance.FSreconditioning','>',0);
			end
			% CHECK THIS ONLY WORKS FOR MOLHO
			if md.flowequation.isMOLHO
				md = checkfield(md,'fieldname','stressbalance.spcvx_base','Inf',1,'timeseries',1);
				md = checkfield(md,'fieldname','stressbalance.spcvy_base','Inf',1,'timeseries',1);
				md = checkfield(md,'fieldname','stressbalance.spcvx_shear','Inf',1,'timeseries',1);
				md = checkfield(md,'fieldname','stressbalance.spcvy_shear','Inf',1,'timeseries',1);
			end
		end % }}}
		function list=defaultoutputs(self,md) % {{{

			if dimension(md.mesh)==3,
				list = {'Vx','Vy','Vz','Vel','Pressure'};
			elseif dimension(md.mesh)==2,
				list = {'Vx','Vy','Vel','Pressure'};
			else
				error('mesh type not supported yet');
			end

		end % }}}
		function disp(self) % {{{

			disp(sprintf('   StressBalance solution parameters:'));

			disp(sprintf('\n      %s','Convergence criteria:'));
			fielddisplay(self,'restol','mechanical equilibrium residual convergence criterion');
			fielddisplay(self,'reltol','velocity relative convergence criterion, NaN: not applied');
			fielddisplay(self,'abstol','velocity absolute convergence criterion, NaN: not applied');
			fielddisplay(self,'isnewton','0: Picard''s fixed point, 1: Newton''s method, 2: hybrid');
			fielddisplay(self,'maxiter','maximum number of nonlinear iterations');

			disp(sprintf('\n      %s','boundary conditions:'));
			fielddisplay(self,'spcvx','x-axis velocity constraint (NaN means no constraint) [m/yr]');
			fielddisplay(self,'spcvy','y-axis velocity constraint (NaN means no constraint) [m/yr]');
			fielddisplay(self,'spcvz','z-axis velocity constraint (NaN means no constraint) [m/yr]');

			disp(sprintf('\n      %s','MOLHO boundary conditions:'));
			fielddisplay(self,'spcvx_base','x-axis basal velocity constraint (NaN means no constraint) [m/yr]');
			fielddisplay(self,'spcvy_base','y-axis basal velocity constraint (NaN means no constraint) [m/yr]');
			fielddisplay(self,'spcvx_shear','x-axis shear velocity constraint (NaN means no constraint) [m/yr]');
			fielddisplay(self,'spcvy_shear','y-axis shear velocity constraint (NaN means no constraint) [m/yr]');

			disp(sprintf('\n      %s','Rift options:'));
			fielddisplay(self,'rift_penalty_threshold','threshold for instability of mechanical constraints');
			fielddisplay(self,'rift_penalty_lock','number of iterations before rift penalties are locked');

			disp(sprintf('\n      %s','Penalty options:'));
			fielddisplay(self,'penalty_factor','offset used by penalties: penalty = Kmax*10^offset');
			fielddisplay(self,'vertex_pairing','pairs of vertices that are penalized');

			disp(sprintf('\n      %s','Other:'));
			fielddisplay(self,'shelf_dampening','use dampening for floating ice ? Only for FS model');
			fielddisplay(self,'FSreconditioning','multiplier for incompressibility equation. Only for FS model');
			fielddisplay(self,'referential','local referential');
			fielddisplay(self,'loadingforce','loading force applied on each point [N/m^3]');
			fielddisplay(self,'requested_outputs','additional outputs requested');

		end % }}}
		function marshall(self,prefix,md,fid) % {{{

			WriteData(fid,prefix,'object',self,'class','stressbalance','fieldname','vertex_pairing','format','DoubleMat','mattype',3);

			yts=md.constants.yts;

			WriteData(fid,prefix,'object',self,'class','stressbalance','fieldname','spcvx','format','DoubleMat','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'class','stressbalance','fieldname','spcvy','format','DoubleMat','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'class','stressbalance','fieldname','spcvz','format','DoubleMat','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',self,'class','stressbalance','fieldname','restol','format','Double');
			WriteData(fid,prefix,'object',self,'class','stressbalance','fieldname','reltol','format','Double');
			WriteData(fid,prefix,'object',self,'class','stressbalance','fieldname','abstol','format','Double','scale',1./yts);
			WriteData(fid,prefix,'object',self,'class','stressbalance','fieldname','isnewton','format','Integer');
			WriteData(fid,prefix,'object',self,'class','stressbalance','fieldname','FSreconditioning','format','Double');
			WriteData(fid,prefix,'object',self,'class','stressbalance','fieldname','maxiter','format','Integer');
			WriteData(fid,prefix,'object',self,'class','stressbalance','fieldname','shelf_dampening','format','Integer');
			WriteData(fid,prefix,'object',self,'class','stressbalance','fieldname','penalty_factor','format','Double');
			WriteData(fid,prefix,'object',self,'class','stressbalance','fieldname','rift_penalty_lock','format','Integer');
			WriteData(fid,prefix,'object',self,'class','stressbalance','fieldname','rift_penalty_threshold','format','Integer');
			WriteData(fid,prefix,'object',self,'class','stressbalance','fieldname','referential','format','DoubleMat','mattype',1);

			if size(self.loadingforce,2)==3,
				WriteData(fid,prefix,'data',self.loadingforce(:,1),'format','DoubleMat','mattype',1,'name','md.stressbalance.loadingforcex');
				WriteData(fid,prefix,'data',self.loadingforce(:,2),'format','DoubleMat','mattype',1,'name','md.stressbalance.loadingforcey');
				WriteData(fid,prefix,'data',self.loadingforce(:,3),'format','DoubleMat','mattype',1,'name','md.stressbalance.loadingforcez');
			end

			%process requested outputs
			outputs = self.requested_outputs;
			pos  = find(ismember(outputs,'default'));
			if ~isempty(pos),
				outputs(pos) = [];                         %remove 'default' from outputs
				outputs      = [outputs defaultoutputs(self,md)]; %add defaults
			end
			WriteData(fid,prefix,'data',outputs,'name','md.stressbalance.requested_outputs','format','StringArray');
			% for MOLHO
			if (md.flowequation.isMOLHO)
				WriteData(fid,prefix,'object',self,'class','stressbalance','fieldname','spcvx_base','format','DoubleMat','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
				WriteData(fid,prefix,'object',self,'class','stressbalance','fieldname','spcvy_base','format','DoubleMat','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
				WriteData(fid,prefix,'object',self,'class','stressbalance','fieldname','spcvx_shear','format','DoubleMat','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
				WriteData(fid,prefix,'object',self,'class','stressbalance','fieldname','spcvy_shear','format','DoubleMat','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			end
		end % }}}
		function savemodeljs(self,fid,modelname) % {{{

			writejs1Darray(fid,[modelname '.stressbalance.spcvx'],self.spcvx);
			writejs1Darray(fid,[modelname '.stressbalance.spcvy'],self.spcvy);
			writejs1Darray(fid,[modelname '.stressbalance.spcvz'],self.spcvz);
			writejsdouble(fid,[modelname '.stressbalance.restol'],self.restol);
			writejsdouble(fid,[modelname '.stressbalance.reltol'],self.reltol);
			writejsdouble(fid,[modelname '.stressbalance.abstol'],self.abstol);
			writejsdouble(fid,[modelname '.stressbalance.isnewton'],self.isnewton);
			writejsdouble(fid,[modelname '.stressbalance.FSreconditioning'],self.FSreconditioning);
			writejsdouble(fid,[modelname '.stressbalance.maxiter'],self.maxiter);
			writejsdouble(fid,[modelname '.stressbalance.shelf_dampening'],self.shelf_dampening);
			writejs1Darray(fid,[modelname '.stressbalance.vertex_pairing'],self.vertex_pairing);
			writejsdouble(fid,[modelname '.stressbalance.penalty_factor'],self.penalty_factor);
			writejsdouble(fid,[modelname '.stressbalance.rift_penalty_lock'],self.rift_penalty_lock);
			writejsdouble(fid,[modelname '.stressbalance.rift_penalty_threshold'],self.rift_penalty_threshold);
			writejs2Darray(fid,[modelname '.stressbalance.referential'],self.referential);
			writejs2Darray(fid,[modelname '.stressbalance.loadingforce'],self.loadingforce);
			writejscellstring(fid,[modelname '.stressbalance.requested_outputs'],self.requested_outputs);

			writejs1Darray(fid,[modelname '.stressbalance.spcvx_base'],self.spcvx_shear);
			writejs1Darray(fid,[modelname '.stressbalance.spcvy_base'],self.spcvy_shear);
			writejs1Darray(fid,[modelname '.stressbalance.spcvx_shear'],self.spcvx_shear);
			writejs1Darray(fid,[modelname '.stressbalance.spcvy_shear'],self.spcvy_shear);
		end % }}}
	end
end
