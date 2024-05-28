//TIMESTEPPINGADAPTIVE class definition
//
//   Usage:
//      timesteppingadaptive=new timesteppingadaptive();

function timesteppingadaptive (){
	//methods
	this.setdefaultparameters = function(){// {{{
		//time between 2 time steps
		this.time_step_min=0.01;
		this.time_step_max=10.;

		//final time
		this.final_time=10.*this.time_step_max;

		//time adaptation? 
		this.cfl_coefficient=0.5;

		//should we interpolate forcings between timesteps?
		this.interp_forcings=1;
		this.cycle_forcing=0;
	}// }}}
	this.disp= function(){// {{{

		var unit;
		console.log(sprintf('   timesteppingadaptive parameters:'));
		unit = 'yr';
		fielddisplay(this,'start_time','simulation starting time ['+ unit + ']');
		fielddisplay(this,'final_time','final time to stop the simulation ['+ unit + ']');
		fielddisplay(this,'time_step_min','minimum length of time steps [' +unit+ ']');
		fielddisplay(this,'time_step_max','maximum length of time steps [' +unit+ ']');
		fielddisplay(this,'cfl_coefficient','coefficient applied to cfl condition');
		fielddisplay(this,'interp_forcings','interpolate in time between requested forcing values ? (0 or 1)');
		fielddisplay(this,'cycle_forcing','cycle through forcing ? (0 or 1)');
		fielddisplay(this,'coupling_time','coupling time steps with ocean model [' +unit+ ']');

	}// }}}
	this.classname= function(){// {{{
		return "timesteppingadaptive";

	}// }}}
	this.checkconsistency = function(md,solution,analyses) { //{{{

		checkfield(md,'fieldname','timestepping.start_time','numel',[1],'NaN',1,'Inf',1);
		checkfield(md,'fieldname','timestepping.final_time','numel',[1],'NaN',1,'Inf',1);
		checkfield(md,'fieldname','timestepping.time_step_min','numel',[1],'>=',0,'NaN',1,'Inf',1);
		checkfield(md,'fieldname','timestepping.time_step_max','numel',[1],'>=',md.timestepping.time_step_max,'NaN',1,'Inf',1);
		checkfield(md,'fieldname','timestepping.cfl_coefficient','numel',[1],'>',0,'<=',1);
		checkfield(md,'fieldname','timestepping.interp_forcings','numel',[1],'values',[0,1]);
		checkfield(md,'fieldname','timestepping.cycle_forcing','numel',[1],'values',[0,1]);
		if (this.final_time-this.start_time<0){
			md.checkmessage('timestepping.final_time should be larger than timestepping.start_time');
		}
		checkfield(md,'fieldname','timestepping.coupling_time','numel',[1],'>=',0,'NaN',1,'Inf',1);
	} // }}}
	this.marshall=function(md,prefix,fid) { //{{{

		var scale;
		scale = md.constants.yts;

		WriteData(fid,prefix,'name','md.timestepping.type','data',2,'format','Integer');
		WriteData(fid,prefix,'object',this,'class','timestepping','fieldname','start_time','format','Double','scale',scale);
		WriteData(fid,prefix,'object',this,'class','timestepping','fieldname','final_time','format','Double','scale',scale);
		WriteData(fid,prefix,'object',this,'class','timestepping','fieldname','time_step_min','format','Double','scale',scale);
		WriteData(fid,prefix,'object',this,'class','timestepping','fieldname','time_step_max','format','Double','scale',scale);
		WriteData(fid,prefix,'object',this,'class','timestepping','fieldname','cfl_coefficient','format','Double');
		WriteData(fid,prefix,'object',this,'class','timestepping','fieldname','interp_forcings','format','Boolean');
		WriteData(fid,prefix,'object',this,'class','timestepping','fieldname','cycle_forcing','format','Boolean');
		WriteData(fid,prefix,'object',this,'class','timestepping','fieldname','coupling_time','format','Double','scale',scale);

	}//}}}
	this.fix=function() { //{{{
	}//}}}
	//properties 
	// {{{
	this.start_time      = 0.;
	this.final_time      = 0.;
	this.time_step_min   = 0.;
	this.time_step_max   = 0.;
	this.cfl_coefficient = 0.;
	this.interp_forcings = 1;
	this.cycle_forcing   = 0;
	this.coupling_time   = 0.;

	this.setdefaultparameters();
	//}}}
}
