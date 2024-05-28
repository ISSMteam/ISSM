//TRANSIENT class definition
//
//   Usage:
//      transient =new transient();

function transient (){
	//methods
	this.setdefaultparameters = function(){// {{{

		//full analysis: Stressbalance, Masstransport and Thermal but no groundingline migration for now
		this.isage             = 0;
		this.issmb             = 1;
		this.ismasstransport   = 1;
		this.isstressbalance   = 1;
		this.isthermal         = 0;
		this.isgroundingline   = 0;
		this.isgia             = 0;
		this.isesa             = 0;
		this.isdamageevolution = 0;
		this.ismovingfront     = 0;
		this.ishydrology       = 0;
		this.isdebris          = 0;
		this.issampling        = 0;
		this.isslr             = 0;
		this.isoceancoupling   = 0;
		this.iscoupler         = 0;
		this.amr_frequency     = 0;

		//default output
		this.requested_outputs=['default'];

	}// }}}
	this.disp= function(){// {{{

		console.log(sprintf('   transient solution parameters:'));

		fielddisplay(this,'isage','indicates whether age model is requested in the transient');
		fielddisplay(this,'issmb','indicates whether a surface mass balance solution is used in the transient');
		fielddisplay(this,'ismasstransport','indicates whether a masstransport solution is used in the transient');
		fielddisplay(this,'isstressbalance','indicates whether a stressbalance solution is used in the transient');
		fielddisplay(this,'isthermal','indicates whether a thermal solution is used in the transient');
		fielddisplay(this,'isgroundingline','indicates whether a groundingline migration is used in the transient');
		fielddisplay(this,'isgia','indicates whether a postglacial rebound model is used in the transient');
		fielddisplay(this,'isesa','indicates whether an elastic adjustment model is used in the transient');
		fielddisplay(this,'isdamageevolution','indicates whether damage evolution is used in the transient');
		fielddisplay(this,'ismovingfront','indicates whether a moving front capability is used in the transient');
		fielddisplay(this,'ishydrology','indicates whether an hydrology model is used');
		fielddisplay(this,'isslr','indicates whether a sea-level rise model is used');
		fielddisplay(this,'isoceancoupling','indicates whether a coupling with an ocean model is used in the transient');
		fielddisplay(this,'iscoupler','indicates whether different models are being run with need for coupling');
		fielddisplay(this,'amr_frequency','frequency at which mesh is refined in simulations with multiple time_steps');
		fielddisplay(this,'requested_outputs','list of additional outputs requested');


	}// }}}
	this.classname= function(){// {{{
		return "transient";
	}// }}}
		this.checkconsistency = function(md,solution,analyses) { // {{{

			//Early return
			if (solution!='TransientSolution') return;

			checkfield(md,'fieldname','transient.isage','numel',[1],'values',[0, 1]);
			checkfield(md,'fieldname','transient.issmb','numel',[1],'values',[0, 1]);
			checkfield(md,'fieldname','transient.ismasstransport','numel',[1],'values',[0, 1]);
			checkfield(md,'fieldname','transient.isstressbalance','numel',[1],'values',[0, 1]);
			checkfield(md,'fieldname','transient.isthermal','numel',[1],'values',[0, 1]);
			checkfield(md,'fieldname','transient.isgroundingline','numel',[1],'values',[0, 1]);
			checkfield(md,'fieldname','transient.isgia','numel',[1],'values',[0, 1]);
			checkfield(md,'fieldname','transient.isesa','numel',[1],'values',[0, 1]);
			checkfield(md,'fieldname','transient.isdamageevolution','numel',[1],'values',[0, 1]);
			checkfield(md,'fieldname','transient.ismovingfront','numel',[1],'values',[0, 1]);
			checkfield(md,'fieldname','transient.ishydrology','numel',[1],'values',[0, 1]);
			checkfield(md,'fieldname','transient.isslr','numel',[1],'values',[0, 1]);
			checkfield(md,'fieldname','transient.isoceancoupling','numel',[1],'values',[0, 1]);
			checkfield(md,'fieldname','transient.iscoupler','numel',[1],'values',[0, 1]);
			checkfield(md,'fieldname','transient.amr_frequency','numel',[1],'>=',0,'NaN',1,'Inf',1);
			checkfield(md,'fieldname','transient.requested_outputs','stringrow',1);
		} // }}}
		this.marshall=function(md,prefix,fid) { //{{{

			prefix='md.transient';
			WriteData(fid,prefix,'object',this,'fieldname','isage','format','Boolean');
			WriteData(fid,prefix,'object',this,'fieldname','issmb','format','Boolean');
			WriteData(fid,prefix,'object',this,'fieldname','ismasstransport','format','Boolean');
			WriteData(fid,prefix,'object',this,'fieldname','isoceantransport','format','Boolean');
			WriteData(fid,prefix,'object',this,'fieldname','isstressbalance','format','Boolean');
			WriteData(fid,prefix,'object',this,'fieldname','isthermal','format','Boolean');
			WriteData(fid,prefix,'object',this,'fieldname','isgroundingline','format','Boolean');
			WriteData(fid,prefix,'object',this,'fieldname','isesa','format','Boolean');
			WriteData(fid,prefix,'object',this,'fieldname','isdamageevolution','format','Boolean');
			WriteData(fid,prefix,'object',this,'fieldname','ishydrology','format','Boolean');
			WriteData(fid,prefix,'object',this,'fieldname','ismovingfront','format','Boolean');
			WriteData(fid,prefix,'object',this,'fieldname','isdebris','format','Boolean');
			WriteData(fid,prefix,'object',this,'fieldname','issampling','format','Boolean');
			WriteData(fid,prefix,'object',this,'fieldname','isslc','format','Boolean');
			WriteData(fid,prefix,'object',this,'fieldname','isoceancoupling','format','Boolean');
			WriteData(fid,prefix,'object',this,'fieldname','amr_frequency','format','Integer');

			//process requested outputs
			var outputs = this.requested_outputs;
			for (var i=0;i<outputs.length;i++){
				if (outputs[i] == 'default') {
					outputs.splice(i,1);
					var newoutputs=this.defaultoutputs(md);
					for (var j=0;j<newoutputs.length;j++) outputs.push(newoutputs[j]);
				}
			}			
			WriteData(fid,prefix,'data',outputs,'name','md.transient.requested_outputs','format','StringArray');
		}//}}}
		this.defaultoutputs = function(md) { //{{{
			return [];
		}//}}}
		this.fix=function() { //{{{
		}//}}}
	//properties 
	// {{{

	this.isage             = 0;
	this.issmb             = 0;
	this.ismasstransport   = 0;
	this.isstressbalance   = 0;
	this.isthermal         = 0;
	this.isgroundingline   = 0;
	this.isgia             = 0;
	this.isesa             = 0;
	this.isdamageevolution = 0;
	this.ismovingfront     = 0;
	this.ishydrology       = 0;
	this.isslr             = 0;
	this.isoceancoupling   = 0;
	this.iscoupler         = 0;
	this.amr_frequency	   = 0;
	this.requested_outputs = [];

	this.setdefaultparameters();
	//}}}
}
