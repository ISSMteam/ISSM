//RIFTS class definition
//
//   Usage:
//      rifts=new rifts();

function rifts (){
	//methods
	this.setdefaultparameters = function(){// {{{
	}// }}}
	this.classname= function(){// {{{
		return "rifts";
	}// }}}
	this.disp= function(){// {{{
		console.log(sprintf('   rifts class echo:'));
		fielddisplay(this,'riftstruct','structure containing all rift information (vertices coordinates, segments, type of melange, ...)');
		fielddisplay(this,'riftproperties','');
	}// }}}
		this.checkconsistency = function(md,solution,analyses) { //{{{
			var numrifts;
			if (isNaN(this.riftstruct) | this.riftstruct.length==0){
				numrifts=0;
			}
			else{
				numrifts=this.riftstruct.length;
			}
			if (numrifts){
				if (!(md.mesh.domaintype() == '2Dhorizontal')){
					md.checkmessage('models with rifts are only supported in 2d for now!');
				}
				if (!IsArray(this.riftstruct)){
					md.checkmessage('rifts.riftstruct should be a structure!');
				}
				for(var i=0;i<md.mesh.segmentmarkers.length;i++){
					if (md.mesh.segmentmarkers[i]>=2){
						//We have segments with rift markers, but no rift structure!
						md.checkmessage(['model should be processed for rifts (run meshprocessrifts)!']);
						break;
					}
				}
				for (var i=0;i<numrifts;i++){
					md = checkfield(md,'fieldname',sprintf('rifts.riftstruct[%i].fill',i),'values',['Water', 'Air', 'Ice', 'Melange']);
				}
			}
			else{
				if (!isNaN(this.riftstruct)) md.checkmessage('riftstruct should be NaN since numrifts is 0!');
			}
		} //}}}
		this.marshall=function(md,prefix,fid) { //{{{

			var numrifts;
			//Process rift info
			if ((this.riftstruct.length==0) | (this.riftstruct == null)){
				numrifts=0;
			}
			else{
				numrifts=this.riftstruct.length;
			}
			var numpairs=0;
			for (var i=0;i<numrifts;i++){
				numpairs=numpairs+this.riftstruct[i].penaltypairs.length;
			}

			// 2 for nodes + 2 for elements+ 2 for  normals + 1 for length + 1 for fill + 1 for friction + 1 for fraction + 1 for fractionincrement + 1 for state.
			data=Create2DArray(numpairs,12);
			var count=0;
			for (var i=0;i<numrifts;i++){
				numpairsforthisrift=this.riftstruct[i].penaltypairs.length;
				for(var j=0;j<numpairsforthisrift;j++){
					for(var k=0;k<7;k++)data[count+j][k]=this.riftstruct[i].penaltypairs;
					data[count+j][7]=this.riftstruct[i].fill;
					data[count+j][8]=this.riftstruct[i].friction;
					data[count+j][9]=this.riftstruct[i].fraction;
					data[count+j][10]=this.riftstruct[i].fractionincrement;
					data[count+j][11]=this.riftstruct[i].state;
					count+=numpairsforthisrift;
				}
			}
			WriteData(fid,prefix,'data',numrifts,'name','md.rifts.numrifts','format','Integer');
			WriteData(fid,prefix,'data',data,'name','md.rifts.riftstruct','format','StringArray');
		}//}}}
		this.fix=function() { //{{{
			this.riftstruct=NullFix(this.riftstruct,'');
			this.riftproperties=NullFix(this.riftproperties,NaN);
		}//}}}
	//properties 
	// {{{
	this.riftstruct     = NaN;
	this.riftproperties = NaN;
	this.setdefaultparameters();
	//}}}
}
