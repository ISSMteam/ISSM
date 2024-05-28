//STEADYSTATE class definition
//
//   Usage:
//      steadystate=new steadystate();

function steadystate (){
	//methods
	this.setdefaultparameters = function(){// {{{

		//maximum of steady state iterations
		this.maxiter=100;

		//Relative tolerance for the steadystate convertgence
		this.reltol=0.01;

		//default output
		this.requested_outputs=['default'];


	}// }}}
	this.disp= function(){// {{{

		console.log(sprintf('   steadystate solution parameters:'));

		fielddisplay(this,'reltol','relative tolerance criterion');
		fielddisplay(this,'maxiter','maximum number of iterations');
		fielddisplay(this,'requested_outputs','additional requested outputs');

	}// }}}
	this.classname= function(){// {{{
		return "steadystate";

	}// }}}
	this.checkconsistency = function(md,solution,analyses) {// {{{

		//Early return
		if (solution!='SteadystateSolution') return;

		if (md.timestepping.time_step!=0){
			md.checkmessage('for a steadystate computation, timestepping.time_step must be zero.');
		}
		checkfield(md,'fieldname','steadystate.requested_outputs','stringrow',1);

		if (isNaN(md.stressbalance.reltol)){
			md.checkmessage('for a steadystate computation, stressbalance.reltol (relative convergence criterion) must be defined!');
		}
	} // }}}
		this.marshall=function(md,prefix,fid) { //{{{
			WriteData(fid,prefix,'object',this,'fieldname','reltol','format','Double');
			WriteData(fid,prefix,'object',this,'fieldname','maxiter','format','Integer');

			//process requested outputs
			var outputs = this.requested_outputs;
			for (var i=0;i<outputs.length;i++){
				if (outputs[i] == 'default') {
					outputs.splice(i,1);
					var newoutputs=this.defaultoutputs(md);
					for (var j=0;j<newoutputs.length;j++) outputs.push(newoutputs[j]);
				}
			}
			WriteData(fid,prefix,'data',outputs,'name','md.steadystate.requested_outputs','format','StringArray');
		}//}}}
		this.defaultoutputs = function(md) { //{{{

			var list=[];

			for (var i=0;i<md.stressbalance.defaultoutputs(md).length;i++)list.push(md.stressbalance.defaultoutputs(md)[i]);
			for (var i=0;i<md.thermal.defaultoutputs(md).length;i++)list.push(md.thermal.defaultoutputs(md)[i]);

			return list;

		}//}}}
		this.fix=function() { //{{{
		}//}}}
	//properties 
	// {{{

	this.reltol            = 0;
	this.maxiter           = 0;
	this.requested_outputs = [];

	this.setdefaultparameters();
	//}}}
}
