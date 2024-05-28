//TIMESTEPPING class definition
//
//   Usage:
//      timestepping=new timestepping();

function timestepping (){
	//methods
	this.setdefaultparameters = function(){ //{{{
		//time between 2 time steps
		this.time_step=1./2.;

		//final time
		this.final_time=10.*this.time_step;

		//should we interpolate forcings between timesteps?
		this.interp_forcing=1;
		this.cycle_forcing=0;
	} //}}}
	this.disp= function(){ //{{{

		var unit;
		console.log(sprintf('   timestepping parameters:'));
		unit = 'yr';
		fielddisplay(this,'start_time','simulation starting time ['+ unit + ']');
		fielddisplay(this,'final_time','final time to stop the simulation ['+ unit + ']');
		fielddisplay(this,'time_step','length of time steps [' +unit+ ']');
		fielddisplay(this,'interp_forcing','interpolate in time between requested forcing values ? (0 or 1)');
		fielddisplay(this,'cycle_forcing','cycle through forcing ? (0 or 1)');
		fielddisplay(this,'coupling_time','length of coupling time steps with ocean model [' +unit+ ']');

	} //}}}
	this.classname= function(){ //{{{
		return "timestepping";

	} //}}}
	this.checkconsistency = function(md,solution,analyses) { //{{{

		checkfield(md,'fieldname','timestepping.start_time','numel',[1],'NaN',1,'Inf',1);
		checkfield(md,'fieldname','timestepping.final_time','numel',[1],'NaN',1,'Inf',1);
		checkfield(md,'fieldname','timestepping.time_step','numel',[1],'>=',0,'NaN',1,'Inf',1);
		checkfield(md,'fieldname','timestepping.interp_forcing','numel',[1],'values',[0,1]);
		checkfield(md,'fieldname','timestepping.cycle_forcing','numel',[1],'values',[0,1]);
		checkfield(md,'fieldname','timestepping.coupling_time','numel',[1],'>=',0,'NaN',1,'Inf',1);
		if (this.final_time-this.start_time<0){
			md.checkmessage('timestepping.final_time should be larger than timestepping.start_time');
		}
		if (solution=='TransientSolution'){
			checkfield(md,'fieldname','timestepping.time_step','numel',[1],'>',0,'NaN',1,'Inf',1);
		}
	} // }}}
	this.marshall=function(md,prefix,fid) { //{{{

		var scale;
		scale = md.constants.yts;

		WriteData(fid,prefix,'name','md.timestepping.type','data',1,'format','Integer');
		WriteData(fid,prefix,'object',this,'fieldname','start_time','format','Double','scale',scale);
		WriteData(fid,prefix,'object',this,'fieldname','final_time','format','Double','scale',scale);
		WriteData(fid,prefix,'object',this,'fieldname','time_step','format','Double','scale',scale);
		WriteData(fid,prefix,'object',this,'fieldname','interp_forcing','format','Boolean');
		WriteData(fid,prefix,'object',this,'fieldname','cycle_forcing','format','Boolean');
		WriteData(fid,prefix,'object',this,'fieldname','coupling_time','format','Double','scale',scale);

	}//}}}
	this.fix=function() { //{{{
	}//}}}
	//properties 
	// {{{
	this.start_time      = 0.;
	this.final_time      = 0.;
	this.time_step       = 0.;
	this.interp_forcing  = 1;
	this.cycle_forcing   = 1;
	this.coupling_time   = 0.;

	this.setdefaultparameters();
	//}}}
}
