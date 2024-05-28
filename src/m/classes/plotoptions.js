//PLOTOPTIONS class definition
//
//   Usage:
//      plotoptions = plotoptions(varargin)

function plotoptions(args) {
	//methods
	this.disp = function (){ // {{{
		console.log(sprintf('\nplotoptions = \n'));
		console.log(sprintf('   figurenumber: %i',this.figurenumber));
		console.log(sprintf('   numberofplots: %i',this.numberofplots));
		if (this.list.length){
			for (var i=0;i<this.list.length;i++){
				console.log(sprintf('\n   options of plot number %i',i+1));
				this.list[i].disp();
			}
		}
		else{
			console.log(sprintf('   list: empty'));
		}
	}
	//}}}
	this.constructor = function (args){ // {{{

		//check length of input
		if (args.length % 2){
			for (i=0;i<args.length;i+=2){
				if (!(typeof args[i] === 'string')){
					console.log('Last valid option: ' + args[i-2]);
					break;
				}
			}
			throw Error('plotoptions error message: invalid parameter/value pair arguments');
		}

		//go through varargin and build list (like pairoptions)
		var rawoptions=new pairoptions(args);
		numoptions=rawoptions.numoptions();

		var counter=0;
		for (i=0;i<numoptions;i++){
			if(typeof args[2*i] === 'string')counter++;
		}
		rawlist=Create2DArray(counter,2);
		var counter=0;
		for (i=0;i<numoptions;i++){
			optionname=args[2*i];
			optionval=args[2*i+1];
			if(typeof optionname === 'string'){
				rawlist[counter][0]=optionname;
				rawlist[counter][1]=optionval;
				counter++;
			}
			else{
				//option is not a string, ignore it
				console.log(sprintf("%s%i%s\n",'WARNING: option number ',i,' is not a string, it will be ignored'));
				rawlist[counter]=[];
				continue
			}
		}
		
			
		//get number of data to be plotted
		numberofplots=rawoptions.fieldoccurrences('data');
		this.numberofplots=numberofplots;

		//figure out wether alloptions flog is on
		if (rawoptions.getfieldvalue('alloptions','off') === 'on') allflag=1;
		else allflag=0;

		//initialize list
		var list=new Array(numberofplots);
		for (i=0;i<numberofplots;i++){
			list[i]=new pairoptions([]);
		}
				
		//process plot options
		for(var i=0;i<rawlist.length;i++){

			//If alloptions flag has is on, apply to all plots
			if (allflag & !(rawlist[i][0] === 'data') & (rawlist[i][0].indexOf('#') == -1)){
				for(var j=0;j<numberofplots;j++){
					list[j].addfield(rawlist[i][0],rawlist[i][1]);
				}
			}
			else if (rawlist[i][0].indexOf('#') != -1){ //option contains '#'

				//get suplot(s) associated
				string=rawlist[i][0].split('#');
				plotnums=string[1];
				field=string[0];

				//divide plotnums if there is a comma ','
				plotnums=plotnums.split(',');

				//loop over plotnums
				for (k=0;k<plotnums.length;k++){
					plotnum=plotnums[k];

					//Empty
					if (plotnum === '') continue;

					else if (plotnum === 'all'){ //pound all
						for(var j=0;j<numberofplots;j++){
							list[j].addfield(field,rawlist[i][1]);
						}
					}
					else if (plotnum.indexOf('-')!=-1){  //pound i-j
						nums=plotnum.split('-');
						if (nums.length!=2) continue;
						if ((nums[0] == '') | (nums[1] === '')){
							throw Error(sprintf("%s%s\n",'the option #i-j is not set properly for ',field));
						}
						for (j=(Number(nums[0])-1);j<(Number(nums[1])); j++){
							list[j].addfield(field,rawlist[i][1]);
						}
					}
					else{ //pound i
						//assign to subplot
						if (Number(plotnum)>numberofplots){
							throw Error(sprintf("%s%s%s%i%s\n",'plotoptions error message: ',field,' cannot be assigned (',plotnum,' exceeds maximum number of plot)'));
						}
						list[Number(plotnum)-1].addfield(field,rawlist[i][1]);
					}
				}
			}
			else{ //assign option field to corresponding subplot

				
				//go through all subplot and assign to the first one free
				var inc=0;
				
				while (inc<numberofplots){
					
					if (!list[inc].exist(rawlist[i][0])){
						list[inc].addfield(rawlist[i][0],rawlist[i][1]);
						break
					}
					else inc++;
				}

				if (inc>numberofplots-1){
					console.log(sprintf("%s%s%s\n",'plot info message: too many ',rawlist[i][0],' options'));
				}
			}
		}

		//check that there is no duplicates
		for (var i=0;i<numberofplots;i++) list[i].deleteduplicates();

		//allocate canvasid automatically
		for (var i=0;i<numberofplots;i++) {
			if (!list[i].exist('canvasid')) {
				list[i].addfield('canvasid',i);
			}
		}

		//Get figure number (should be in options for subplot 1)
		this.figurenumber=list[0].getfieldvalue('figure',1);
		list[0].removefield('figure',0);

		//asign output
		this.list=list;

	} //}}}
	//properties
	// {{{
	this.numberofplots = 0;
	this.figurenumber  = 1;
	this.list          = [];
	this.constructor(args);
	//}}}
}
