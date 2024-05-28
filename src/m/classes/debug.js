//DEBUG class definition
//
//   Usage:
//      debug=new debug();

function debug (){
	//methods
	this.setdefaultparameters = function(){// {{{
	}// }}}
	this.classname= function(){// {{{
		return "debug";
	}// }}}
	this.disp= function(){// {{{
		console.log(sprintf('   debug parameters:'));
		console.log(sprintf('   debug parameters:'));

		fielddisplay(this,'valgrind','use Valgrind to debug (0 or 1)');
		fielddisplay(this,'gprof','use gnu-profiler to find out where the time is spent');
		fielddisplay(this,'profiling','enables profiling (memory, flops, time)');

	}// }}}
		this.marshall=function(md,prefix,fid) { //{{{
			WriteData(fid,prefix,'object',this,'fieldname','profiling','format','Boolean');
		}//}}}
		this.fix=function() { //{{{
		}//}}}

	//properties 
	// {{{
	this.valgrind = false;
	this.gprof    = false;
	this.profiling = false;
	this.setdefaultparameters();
	//}}}
}
