//RADAROVERLAY class definition
//
//   Usage:
//      radaroverlay=new radaroverlay();

function radaroverlay (){
	//methods
	this.setdefaultparameters = function(){// {{{
	}// }}}
	this.disp= function(){// {{{
		console.log(sprintf('   radaroverlay parameters:'));

		fielddisplay(this,'xlim','corresponding x boundaries[m]');
		fielddisplay(this,'ylim','corresponding y boundaries [m]');
		fielddisplay(this,'outerindex','outer triangulation between mesh and bounding box');
		fielddisplay(this,'outerx','outer triangulation x coordinate between mesh and bounding box');
		fielddisplay(this,'outery','outer triangulation y coordinate between mesh and bounding box');

	}// }}}
	//properties 
	// {{{
	this.xlim   = NaN;
	this.ylim   = NaN;
	this.outerindex   = NaN;
	this.outerx   = NaN;
	this.outery   = NaN;
	this.setdefaultparameters();
	//}}}
}
