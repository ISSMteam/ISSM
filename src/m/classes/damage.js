//DAMAGE class definition
//
//   Usage:
//      damage=new damage();

function damage (){
	//methods
	this.setdefaultparameters = function(){// {{{
		
		//damage parameters: 
		this.isdamage=0;
		this.D=0;
		this.law=0;

		this.max_damage=1-1e-5; //if damage reaches 1, solve becomes singular, as viscosity becomes nil

		//Type of stabilization used
		this.stabilization=4;

		//Maximum number of iterations
		this.maxiter=100;

		//finite element interpolation
		this.elementinterp='P1';

		//damage evolution parameters 
		this.stress_threshold=1.3e5;
		this.kappa=2.8;
		this.healing=0;
		this.c1=0;
		this.c2=0;
		this.c3=0;
		this.c4=0;
		this.equiv_stress=0;

		//output default:
		this.requested_outputs=['default'];

	}// }}}
	this.disp= function(){// {{{
		console.log(sprintf('   Damage:\n'));

		fielddisplay(this,'isdamage','is damage mechanics being used? {true,false}');
		if (this.isdamage){
			fielddisplay(this,'law',"damage law ['0: analytical','1: pralong']");
			fielddisplay(this,'D','damage tensor (scalar)');
			fielddisplay(this,'spcdamage','damage constraints (NaN means no constraint)');
			fielddisplay(this,'max_damage','maximum possible damage (0<=max_damage<1)');

			fielddisplay(this,'stabilization','0: no stabilization, 1: artificial diffusion, 2: SUPG (not working), 4: flux corrected transport');
			fielddisplay(this,'maxiter','maximum number of non linear iterations');
			fielddisplay(this,'elementinterp',"interpolation scheme for finite elements {'P1','P2'}");
			fielddisplay(this,'stress_threshold','stress threshold for damage initiation [Pa]');
			fielddisplay(this,'kappa','ductility parameter for stress softening and damage');
			fielddisplay(this,'c1','damage parameter 1');
			fielddisplay(this,'c2','damage parameter 2');
			fielddisplay(this,'c3','damage parameter 3');
			fielddisplay(this,'c4','damage parameter 4');
			fielddisplay(this,'healing','damage healing parameter');
			fielddisplay(this,'equiv_stress','0: von Mises, 1: max principal');
			fielddisplay(this,'requested_outputs','additional outputs requested');
		}
	}// }}}
    this.extrude = function(md) {//{{{
        this.D=project3d(md,'vector',this.D,'type','node');
        this.spcdamage=project3d(md,'vector',this.spcdamage,'type','node');
        return this;
    }//}}}
	this.classname= function(){// {{{
		return "damage";
	}// }}}
		this.checkconsistency = function(md,solution,analyses) { //{{{
			
			checkfield(md,'fieldname','damage.isdamage','values',[1,0]);
			if (this.isdamage){
				checkfield(md,'fieldname','damage.law','numel',[1],'values',[0,1,2]);
				checkfield(md,'fieldname','damage.D','>=',0,'<=',this.max_damage,'size',[md.mesh.numberofvertices ,1]);
				checkfield(md,'fieldname','damage.spcdamage','Inf',1,'timeseries',1);
				checkfield(md,'fieldname','damage.max_damage','<',1,'>=',0);
				checkfield(md,'fieldname','damage.stabilization','numel',[1],'values',[0, 1, 2, 4]);
				checkfield(md,'fieldname','damage.maxiter','>=0',0);
				checkfield(md,'fieldname','damage.elementinterp','values',['P1','P2']);
				checkfield(md,'fieldname','damage.stress_threshold','>=',0);
				checkfield(md,'fieldname','damage.kappa','>',1);
				checkfield(md,'fieldname','damage.healing','>=',0);
				checkfield(md,'fieldname','damage.c1','>=',0);
				checkfield(md,'fieldname','damage.c2','>=',0);
				checkfield(md,'fieldname','damage.c3','>=',0);
				checkfield(md,'fieldname','damage.c4','>=',0);
				checkfield(md,'fieldname','damage.equiv_stress','numel',[1],'values',[0, 1]);
				checkfield(md,'fieldname','damage.requested_outputs','stringrow',1);
			}
			else if (this.law!=0){
				if (solution=='DamageEvolutionSolution'){
					throw Error('Invalid evolution law (md.damage.law) for a damage solution');
				}
			}
		} //}}}
		this.marshall=function(md,prefix,fid) { //{{{
		
			WriteData(fid,prefix,'object',this,'fieldname','isdamage','format','Boolean');
			if (this.isdamage){
				WriteData(fid,prefix,'object',this,'fieldname','law','format','Integer');
				WriteData(fid,prefix,'object',this,'fieldname','D','format','DoubleMat','mattype',1);
				WriteData(fid,prefix,'object',this,'fieldname','spcdamage','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
				WriteData(fid,prefix,'object',this,'fieldname','max_damage','format','Double');

				WriteData(fid,prefix,'object',this,'fieldname','stabilization','format','Integer');
				WriteData(fid,prefix,'object',this,'fieldname','maxiter','format','Integer');
				WriteData(fid,prefix,'name','md.damage.elementinterp','data',this.elementinterp,'format','String');
				WriteData(fid,prefix,'object',this,'fieldname','stress_threshold','format','Double');
				WriteData(fid,prefix,'object',this,'fieldname','kappa','format','Double');
				WriteData(fid,prefix,'object',this,'fieldname','c1','format','Double');
				WriteData(fid,prefix,'object',this,'fieldname','c2','format','Double');
				WriteData(fid,prefix,'object',this,'fieldname','c3','format','Double');
				WriteData(fid,prefix,'object',this,'fieldname','c4','format','Double');
				WriteData(fid,prefix,'object',this,'fieldname','healing','format','Double');
				WriteData(fid,prefix,'object',this,'fieldname','equiv_stress','format','Integer');
			}

			//process requested outputs
			var outputs = this.requested_outputs;
			for (var i=0;i<outputs.length;i++){
				if (outputs[i] == 'default') {
					outputs.splice(i,1);
					outputs.push(this.defaultoutputs(md));
				}
			}
			if (this.isdamage){
				WriteData(fid,prefix,'data',outputs,'name','md.damage.requested_outputs','format','StringArray');
			}

		}//}}}
		this.fix=function() { //{{{
		}//}}}
		this.defaultoutputs = function(md){ //{{{

			if (md.mesh.domaintype() == '2Dhorizontal') return 'DamageDbar';
			else return 'DamageD';

		}//}}}
	//properties 
	// {{{
	this.isdamage            = 0;
	this.D                   = NaN;
	this.law                 = 0;
	this.spcdamage           = NaN; 
	this.max_damage          = 0;

	//numerical
	this.stabilization       = 0;
	this.maxiter             = 0;
	this.elementinterp       = '';

	//general parameters for evolution law: 
	this.stress_threshold    = 0;
	this.kappa               = 0;
	this.c1                  = 0;
	this.c2                  = 0;
	this.c3                  = 0;
	this.c4                  = 0;
	this.healing             = 0;
	this.equiv_stress		  = 0;
	this.requested_outputs   = [];

	this.setdefaultparameters();
	//}}}
}
