%FRICTIONNN Native 2-32-32-1 basal drag network for horizontal SSA.
classdef frictionnn
	properties
		C = NaN;
		nn_weights = NaN(1,1185);
		nn_normalization = [0 0 1 1 0 1 1e-12];
	end
	methods
		function md = checkconsistency(self,md,solution,analyses)
			if ~ismember('StressbalanceAnalysis',analyses), return; end
			if ~strcmp(md.mesh.domaintype(),'2Dhorizontal') || any(md.flowequation.element_equation~=2)
				md = checkmessage(md,'frictionnn requires a horizontal SSA mesh.');
			end
			if md.inversion.iscontrol && (~md.autodiff.isautodiff || ~isa(md.inversion,'adm1qn3inversion'))
				md = checkmessage(md,'frictionnn inversion requires autodiff and adm1qn3inversion.');
			end
			md = checkfield(md,'fieldname','friction.C','size',[md.mesh.numberofvertices 1],'NaN',1,'Inf',1,'>=',0);
			md = checkfield(md,'fieldname','friction.nn_weights','size',[1 1185],'NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','friction.nn_normalization','size',[1 7],'NaN',1,'Inf',1);
			if any(self.nn_normalization([3 4 6 7])<=0)
				md = checkmessage(md,'frictionnn scales and speed floor must be positive.');
			end
		end
		function disp(self)
			disp('Native neural basal friction law (2-32-32-1 ReLU)');
			fielddisplay(self,'C','nodal friction coefficient [SI]');
			fielddisplay(self,'nn_weights','row of per-neuron weights followed by bias, layer by layer');
			fielddisplay(self,'nn_normalization','[two input centers, two input scales, output center, output scale, speed floor (m/s)]');
		end
		function marshall(self,prefix,md,fid)
			WriteData(fid,prefix,'name','md.friction.law','data',21,'format','Integer');
			WriteData(fid,prefix,'object',self,'fieldname','C','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',self,'fieldname','nn_weights','format','DoubleMat','mattype',3);
			WriteData(fid,prefix,'object',self,'fieldname','nn_normalization','format','DoubleMat','mattype',3);
		end
	end
end
