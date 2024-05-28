//STRESSBALANCE class definition
//
//   Usage:
//      stressbalance=new stressbalance();

function stressbalance (){
	//methods
	this.setdefaultparameters = function(){// {{{

		//maximum of non-linear iterations.
		this.maxiter=100;

		//Convergence criterion: absolute, relative and residual
		this.restol=Math.pow(10,-4); 
		this.reltol=0.01;
		this.abstol=10;

		this.FSreconditioning=Math.pow(10,13);
		this.shelf_dampening=0;

		//Penalty factor applied kappa=max(stiffness matrix)*10^penalty_factor
		this.penalty_factor=3;

		//Stop the iterations of rift if below a threshold
		this.rift_penalty_threshold=0;

		//in some solutions, it might be needed to stop a run when only
		//a few constraints remain unstable. For thermal computation, this
		//parameter is often used.
		this.rift_penalty_lock=10;

		//output default:
		this.requested_outputs=['default'];

	}// }}}
	this.disp= function(){// {{{
		console.log(sprintf('   StressBalance solution parameters:'));

		console.log(sprintf('\n      %s','Convergence criteria:'));
		fielddisplay(this,'restol','mechanical equilibrium residual convergence criterion');
		fielddisplay(this,'reltol','velocity relative convergence criterion, NaN: not applied');
		fielddisplay(this,'abstol','velocity absolute convergence criterion, NaN: not applied');
		fielddisplay(this,'isnewton',"0: Picard's fixed point, 1: Newton's method, 2: hybrid");
		fielddisplay(this,'maxiter','maximum number of nonlinear iterations');

		console.log(sprintf('\n      %s','boundary conditions:'));
		fielddisplay(this,'spcvx','x-axis velocity constraint (NaN means no constraint) [m/yr]');
		fielddisplay(this,'spcvy','y-axis velocity constraint (NaN means no constraint) [m/yr]');
		fielddisplay(this,'spcvz','z-axis velocity constraint (NaN means no constraint) [m/yr]');

		console.log(sprintf('\n      %s','Rift options:'));
		fielddisplay(this,'rift_penalty_threshold','threshold for instability of mechanical constraints');
		fielddisplay(this,'rift_penalty_lock','number of iterations before rift penalties are locked');

		console.log(sprintf('\n      %s','Penalty options:'));
		fielddisplay(this,'penalty_factor','offset used by penalties: penalty = Kmax*10^offset');
		fielddisplay(this,'vertex_pairing','pairs of vertices that are penalized');

		console.log(sprintf('\n      %s','Other:'));
		fielddisplay(this,'shelf_dampening','use dampening for floating ice ? Only for FS model');
		fielddisplay(this,'FSreconditioning','multiplier for incompressibility equation. Only for FS model');
		fielddisplay(this,'referential','local referential');
		fielddisplay(this,'loadingforce','loading force applied on each point [N/m^3]');
		fielddisplay(this,'requested_outputs','additional outputs requested');



	}// }}}
	this.classname= function(){// {{{
		return "stressbalance";
	}// }}}
	this.extrude = function(md) {//{{{
		this.spcvx=project3d(md,'vector',this.spcvx,'type','node');
		this.spcvy=project3d(md,'vector',this.spcvy,'type','node');
		this.spcvz=project3d(md,'vector',this.spcvz,'type','node');
		this.referential=project3d(md,'vector',this.referential,'type','node');
		this.loadingforce=project3d(md,'vector',this.loadingforce,'type','node');
		return this;
	}//}}}
	this.checkconsistency = function(md,solution,analyses) { //{{{

		//Early return
		if(ArrayAnyEqual(ArrayIsMember('StressbalanceAnalysis',analyses),0))return;

		checkfield(md,'fieldname','stressbalance.spcvx','Inf',1,'timeseries',1);
		checkfield(md,'fieldname','stressbalance.spcvy','Inf',1,'timeseries',1);
		checkfield(md,'fieldname','stressbalance.spcvz','Inf',1,'timeseries',1);
		checkfield(md,'fieldname','stressbalance.restol','size',[1, 1],'>',0,'NaN',1,'Inf',1);
		checkfield(md,'fieldname','stressbalance.reltol','size',[1, 1]);
		checkfield(md,'fieldname','stressbalance.abstol','size',[1, 1]);
		checkfield(md,'fieldname','stressbalance.isnewton','numel',[1],'values',[0, 1, 2]);
		checkfield(md,'fieldname','stressbalance.FSreconditioning','size',[1, 1],'NaN',1,'Inf',1);
		checkfield(md,'fieldname','stressbalance.maxiter','size',[1, 1],'>=',1);
		checkfield(md,'fieldname','stressbalance.referential','size',[md.mesh.numberofvertices, 6]);
		checkfield(md,'fieldname','stressbalance.loadingforce','size',[md.mesh.numberofvertices, 3]);
		checkfield(md,'fieldname','stressbalance.requested_outputs','stringrow',1);

		//singular solution
		if(!ArrayAnyNaN(md.stressbalance.spcvx) | !ArrayAnyNaN(md.stressbalance.spcvy) |  !ArrayAnyAboveStrict(md.mask.ocean_levelset,0)){
			md = checkmessage(md,'model is not well posed (singular). You need at least one node with fixed velocity!');
			console.log(sprintf('\n !!! Warning: no spc applied, model might not be well posed if no basal friction is applied, check for solution crash\n'));
		}
		//CHECK THAT EACH LINES CONTAINS ONLY NAN VALUES OR NO NAN VALUES
		for(var i=0;i<md.stressbalance.referential.length;i++){
			var sum=0;
			for(j=0;j<md.stressbalance.referential[0].length;j++)sum+=md.stressbalance.referential[i][j];
			if (sum!=0 & sum!=6){
				md = checkmessage(md,'Each line of stressbalance.referential should contain either only NaN values or no NaN values');
				break;
			}
		}
		//CHECK THAT THE TWO VECTORS PROVIDED ARE ORTHOGONAL
		for(var i=0;i<md.stressbalance.referential.length;i++){
			var sum=0;
			for(j=0;j<md.stressbalance.referential[0].length;j++)sum+=md.stressbalance.referential[i][j];
			if(sum==0){
				var dot=0;
				for(j=0;j<3;j++)dot+=md.stressbalance.referential[i][j]*md.stressbalance.referential[i][j+3];
				dot=Math.abs(dot);
				if(dot>Math.pow(10,-18)){
					md.checkmessage('Vectors in stressbalance.referential (columns 1 to 3 and 4 to 6) must be orthogonal');
					break;
				}
			}
		}
		//CHECK THAT NO rotation specified for FS Grounded ice at base
		if (md.mesh.domaintype() == '3D' & md.flowequation.isFS){
			for(var i=0;i<md.mask.ocean_levelset.length;i++){
				if(md.mask.ocean_levelset[i]>0 & md.mesh.vertexonbase[i]){
					if(!ArrayIsNan(md.stressbalance.referential[i])){
						md.checkmessage('no referential should be specified for basal vertices of grounded ice');
						break;
					}
				}
			}
			checkfield(md,'fieldname','stressbalance.FSreconditioning','>',0);
		}
	} // }}}
	this.marshall=function(md,prefix,fid) { //{{{

		WriteData(fid,prefix,'object',this,'class','stressbalance','fieldname','vertex_pairing','format','DoubleMat','mattype',3);

		var yts=md.constants.yts;

		WriteData(fid,prefix,'object',this,'class','stressbalance','fieldname','spcvx','format','DoubleMat','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
		WriteData(fid,prefix,'object',this,'class','stressbalance','fieldname','spcvy','format','DoubleMat','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
		WriteData(fid,prefix,'object',this,'class','stressbalance','fieldname','spcvz','format','DoubleMat','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
		WriteData(fid,prefix,'object',this,'class','stressbalance','fieldname','restol','format','Double');
		WriteData(fid,prefix,'object',this,'class','stressbalance','fieldname','reltol','format','Double');
		WriteData(fid,prefix,'object',this,'class','stressbalance','fieldname','abstol','format','Double','scale',1./yts);
		WriteData(fid,prefix,'object',this,'class','stressbalance','fieldname','isnewton','format','Integer');
		WriteData(fid,prefix,'object',this,'class','stressbalance','fieldname','FSreconditioning','format','Double');
		WriteData(fid,prefix,'object',this,'class','stressbalance','fieldname','maxiter','format','Integer');
		WriteData(fid,prefix,'object',this,'class','stressbalance','fieldname','shelf_dampening','format','Integer');
		WriteData(fid,prefix,'object',this,'class','stressbalance','fieldname','penalty_factor','format','Double');
		WriteData(fid,prefix,'object',this,'class','stressbalance','fieldname','rift_penalty_lock','format','Integer');
		WriteData(fid,prefix,'object',this,'class','stressbalance','fieldname','rift_penalty_threshold','format','Integer');
		WriteData(fid,prefix,'object',this,'class','stressbalance','fieldname','referential','format','DoubleMat','mattype',1);

		if (size(this.loadingforce,1)==3) {
			WriteData(fid,prefix,'data',ArrayCol(this.loadingforce,0),'format','DoubleMat','mattype',1,'name','md.stressbalance.loadingforcex');
			WriteData(fid,prefix,'data',ArrayCol(this.loadingforce,1),'format','DoubleMat','mattype',1,'name','md.stressbalance.loadingforcey');
			WriteData(fid,prefix,'data',ArrayCol(this.loadingforce,2),'format','DoubleMat','mattype',1,'name','md.stressbalance.loadingforcez');
		}

		//process requested outputs
		var outputs = this.requested_outputs;
		for (var i=0;i<outputs.length;i++){
			if (outputs[i] == 'default') {
				outputs.splice(i,1);
				var newoutputs=this.defaultoutputs(md);
				for (var j=0;j<newoutputs.length;j++) outputs.push(newoutputs[j]);
			}
		}
		WriteData(fid,prefix,'data',outputs,'name','md.stressbalance.requested_outputs','format','StringArray');
	}//}}}
	this.defaultoutputs = function(md){ // {{{

		var list;
		if (md.mesh.dimension() == 3) list = ['Vx','Vy','Vz','Vel','Pressure'];
		else if (md.mesh.dimension()==2) list = ['Vx','Vy','Vel','Pressure'];
		else throw Error('mesh type not supported yet');
		return list;

	}//}}}
	this.fix=function() { //{{{
		this.abstol=NullFix(this.abstol,NaN);
		this.rift_penalty_lock=NullFix(this.rift_penalty_lock,NaN);
		this.referential=NullFix(this.referential,NaN);
		this.loadingforce=NullFix(this.loadingforce,NaN);
		this.spcvx=NullFix(this.spcvx,NaN);
		this.spcvy=NullFix(this.spcvy,NaN);
		this.spcvz=NullFix(this.spcvz,NaN);
		if(this.vertex_pairing=[])this.vertex_pairing=NaN;
	}//}}}
	//properties 
	// {{{
	this.spcvx                  = NaN;
	this.spcvy                  = NaN;
	this.spcvz                  = NaN;
	this.restol                 = 0;
	this.reltol                 = 0;
	this.abstol                 = 0;
	this.isnewton               = 0;
	this.FSreconditioning       = 0;
	this.maxiter                = 0;
	this.shelf_dampening        = 0;
	this.vertex_pairing         = NaN;
	this.penalty_factor         = NaN;
	this.rift_penalty_lock      = NaN;
	this.rift_penalty_threshold = 0;
	this.referential            = NaN;
	this.loadingforce           = NaN;
	this.requested_outputs      = []

	this.setdefaultparameters();
	//}}}
}
