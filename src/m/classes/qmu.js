//QMU class definition
//
//   Usage:
//      qmu=new qmu();

function qmu (){
	//methods
	this.setdefaultparameters = function(){// {{{
	}// }}}
	this.disp= function(){// {{{

		console.log(sprintf('   qmu parameters:'));

		fielddisplay(this,'isdakota','is qmu analysis activated?');
		fielddisplay(this,'output','are we outputting ISSM results, default is 0');
		/*
		for (var i=0;i<this.variables.length;i++){
			console.log(sprintf('         variables%s:  (arrays of each variable class)',...
						string_dim(this.variables,i)));
		}
		fnames=fieldnames(this.variables(i));
		maxlen=0;
		for j=1:numel(fnames)
			maxlen=max(maxlen,length(fnames{j}));
		end

			for j=1:numel(fnames)
				console.log(sprintf(['            %-' num2str(maxlen+1) 's:    [%ix%i]    ''%s'''],...
							fnames{j},size(this.variables.(fnames{j})),class(this.variables.(fnames{j}))));
		end
			end
			for i=1:numel(this.responses)
				console.log(sprintf('         responses%s:  (arrays of each response class)',...
							string_dim(this.responses,i)));
		fnames=fieldnames(this.responses(i));
		maxlen=0;
		for j=1:numel(fnames)
			maxlen=max(maxlen,length(fnames{j}));
		end

			for j=1:numel(fnames)
				console.log(sprintf(['            %-' num2str(maxlen+1) 's:    [%ix%i]    ''%s'''],...
							fnames{j},size(this.responses.(fnames{j})),class(this.responses.(fnames{j}))));
		end
			end
			fielddisplay(this,'numberofresponses','number of responses') 
			for i=1:numel(this.method);
		if strcmp(class(this.method(i)),'dakota_method')
			console.log(sprintf('            method%s :    ''%s''',...
						string_dim(this.method,i),this.method(i).method));
		end
			end
			for i=1:numel(this.params)
				console.log(sprintf('         params%s:  (array of method-independent parameters)',...
							string_dim(this.params,i)));
		fnames=fieldnames(this.params(i));
		maxlen=0;
		for j=1:numel(fnames)
			maxlen=max(maxlen,length(fnames{j}));
		end

			for j=1:numel(fnames)
				console.log(sprintf(['            %-' num2str(maxlen+1) 's: %s'],...
							fnames{j},any2str(this.params(i).(fnames{j}))));
		end
			end
			for i=1:numel(this.results)
				console.log(sprintf('         results%s:  (information from dakota files)',...
							string_dim(this.results,i)));
		fnames=fieldnames(this.results(i));
		maxlen=0;
		for j=1:numel(fnames)
			maxlen=max(maxlen,length(fnames{j}));
		end

			for j=1:numel(fnames)
				console.log(sprintf(['            %-' num2str(maxlen+1) 's:    [%ix%i]    ''%s'''],...
							fnames{j},size(this.results.(fnames{j})),class(this.results.(fnames{j}))));
		end
			end
			fielddisplay(this,'partition','user provided mesh partitioning, defaults to metis if not specified') 
			fielddisplay(this,'numberofpartitions','number of partitions for semi-discrete qmu') 
			fielddisplay(this,'variabledescriptors','');
		fielddisplay(this,'responsedescriptors','');
		fielddisplay(this,'method','array of dakota_method class');
		fielddisplay(this,'mass_flux_profile_directory','directory for mass flux profiles');
		fielddisplay(this,'mass_flux_profiles','list of mass_flux profiles');
		fielddisplay(this,'mass_flux_segments','');
		fielddisplay(this,'adjacency','');
		fielddisplay(this,'vertex_weight','weight applied to each mesh vertex');
		*/

	}// }}}
    this.extrude = function(md) {//{{{
	if (!isNaN(this.partition)) this.partition=project3d(md,'vector',ArrayTranspose(this.partition),'type','node');
        return this;
    }//}}}
	this.classname= function(){// {{{
		return "qmu";
	}// }}}
		this.checkconsistency = function(md,solution,analyses) { //{{{

			///Early return
			if (!md.qmu.isdakota) return;
			else md.checkmessage('qmu runs not supported yet!');

		} // }}}
		this.marshall=function(md,prefix,fid) { //{{{
			WriteData(fid,prefix,'object',this,'fieldname','isdakota','format','Boolean');
			WriteData(fid,prefix,'object',this,'fieldname','output','format','Boolean');
			if (!this.isdakota){
				WriteData(fid,prefix,'data',0,'name','md.qmu.mass_flux_segments_present','format','Boolean');
			}
			else{
				WriteData(fid,prefix,'object',this,'fieldname','partition','format','DoubleMat','mattype',2);
				WriteData(fid,prefix,'object',this,'fieldname','numberofpartitions','format','Integer');
				WriteData(fid,prefix,'object',this,'fieldname','numberofresponses','format','Integer');
				WriteData(fid,prefix,'object',this,'fieldname','variabledescriptors','format','StringArray');
				WriteData(fid,prefix,'object',this,'fieldname','responsedescriptors','format','StringArray');
				if (this.mass_flux_segments.length){
					WriteData(fid,prefix,'data',this.mass_flux_segments,'name','md.qmu.mass_flux_segments','format','MatArray');
					flag=true; 
				}
				else flag=false; 
				WriteData(fid,prefix,'data',flag,'name','md.qmu.mass_flux_segments_present','format','Boolean');
			}
		}//}}}
		this.fix=function() { //{{{
		}//}}}
	//properties 
	// {{{

	this.isdakota                    = 0;
	this.output                      = 0;
	this.variables                   = []
	this.responses                   = [];
	this.method                      = []
	this.params                      = []
	this.results                     = []
	this.partition                   = NaN;
	this.numberofpartitions          = 0;
	this.numberofresponses           = 0;
	this.variabledescriptors         = []
	this.responsedescriptors         = []
	this.mass_flux_profile_directory = NaN;
	this.mass_flux_profiles          = NaN;
	this.mass_flux_segments          = []
	this.adjacency                   = NaN;
	this.vertex_weight               = NaN;

	this.setdefaultparameters();
	//}}}
}
