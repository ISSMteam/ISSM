//ESA class definition
//
//   Usage:
//      esa=esa();

function esa(){
	//methods
		this.setdefaultparameters = function (){ //{{{
		
		//numerical discretization accuracy
		this.degacc=.01;
	
		//computational flags:
		this.hemisphere=0;

		//output default:
		this.requested_outputs=['default'];

		//transitions should be a cell array of vectors: 
		this.transitions=[];
		
		}// }}}
		this.checkconsistency = function(md,solution,analyses) { //{{{

			//Early return
			if(ArrayAnyEqual(ArrayIsMember('EsaAnalysis',analyses),0))return;
			
			md = checkfield(md,'fieldname','esa.deltathickness','NaN',1,'Inf',1,'size',[md.mesh.numberofelements, 1]);
			md = checkfield(md,'fieldname','esa.love_h','NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','esa.love_l','NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','esa.hemisphere','NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','esa.degacc','size',[1, 1],'>=',1e-10);
			md = checkfield(md,'fieldname','esa.requested_outputs','stringrow',1);
			
			//check that love numbers are provided at the same level of accuracy: 
			if (this.love_h.length != this.love_l.length){
				throw Error('esa error message: love numbers should be provided at the same level of accuracy');
			}

		} // }}}
		this.defaultoutputs = function(md){ // {{{
			return ['EsaUmotion'];
		}//}}}
		this.classname= function(){// {{{
			return "esa";
		}// }}}
		this.disp= function(){// {{{
			
		console.log(sprintf('   esa solution parameters:'));

		fielddisplay(this,'deltathickness','thickness change: ice height equivalent [m]');
		fielddisplay(this,'love_h','load Love number for radial displacement');
		fielddisplay(this,'love_l','load Love number for horizontal displacements'); 
		fielddisplay(this,'hemisphere','North-south, East-west components of 2-D horiz displacement vector: -1 south, 1 north'); 
		fielddisplay(this,'degacc',"accuracy (default .01 deg) for numerical discretization of the Green's functions");
		fielddisplay(this,'transitions','indices into parts of the mesh that will be icecaps');
		fielddisplay(this,'requested_outputs','additional outputs requested (default: EsaUmotion)');
		} //}}}
		this.marshall=function(md,prefix,fid) { //{{{

			WriteData(fid,prefix,'object',this,'fieldname','deltathickness','format','DoubleMat','mattype',2);
			WriteData(fid,prefix,'object',this,'fieldname','love_h','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',this,'fieldname','love_l','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',this,'fieldname','hemisphere','format','Integer');
			WriteData(fid,prefix,'object',this,'fieldname','degacc','format','Double');
			WriteData(fid,prefix,'object',this,'fieldname','transitions','format','MatArray');

			//process requested outputs
			var outputs = this.requested_outputs;
			for (var i=0;i<outputs.length;i++){
				if (outputs[i] == 'default') {
					outputs.splice(i,1);
					var newoutputs=this.defaultoutputs(md);
					for (var j=0;j<newoutputs.length;j++) outputs.push(newoutputs[j]);
				}
			}
			WriteData(fid,prefix,'data',outputs,'name','md.esa.requested_outputs','format','StringArray');
		}//}}}
		this.fix=function() { //{{{
			this.deltathickness=NullFix(this.deltathickness,NaN);
			this.love_h=NullFix(this.love_h,NaN);
			this.love_l=NullFix(this.love_l,NaN);
			this.hemisphere=NullFix(this.hemisphere,NaN);
			this.degacc=NullFix(this.degacc,NaN);
		}//}}}
	//properties
	//{{{
	this.deltathickness = NaN;
	this.love_h         = 0; //provided by PREM model
	this.love_l         = 0; //idam
	this.hemisphere     = 0;
	this.degacc         = 0;
	this.requested_outputs = [];
	this.transitions    = [];
	this.setdefaultparameters();
	//}}}
}
