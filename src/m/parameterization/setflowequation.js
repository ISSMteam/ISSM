function setflowequation(md){
//SETFLOWEQUATION - associate a solution type to each element
//
//   This routine works like plotmodel: it works with an even number of inputs
//   'SIA','SSA','L1L2','MOLHO','HO','FS' and 'fill' are the possible options
//   that must be followed by the corresponding exp file or flags list
//   It can either be a domain file (argus type, .exp extension), or an array of element flags. 
//   If user wants every element outside the domain to be 
//   setflowequationd, add '~' to the name of the domain file (ex: '~HO.exp');
//   an empty string '' will be considered as an empty domain
//   a string 'all' will be considered as the entire domain
//   You can specify the type of coupling, 'penalties' or 'tiling', to use with the input 'coupling'
//   NB: L1L2 and MOLHO cannot currently be coupled to any other ice flow model
//
//   Usage:
//      setflowequation(md,varargin)
//
//   Example:
//      setflowequation(md,'HO',HO,'fill','SIA','coupling','tiling');

	//some checks on list of arguments
	if(arguments.length<3) throw Error('setflowequation error message');

	//Process options
	var args = Array.prototype.slice.call(arguments);
	var options = new pairoptions(args.slice(1,args.length));
	options.deleteduplicates(1);

	//Find_out what kind of coupling to use
	coupling_method=options.getfieldvalue('coupling','tiling');
	if ((coupling_method != 'tiling') & !(coupling_method != 'penalties')){
		throw error('coupling type can only be: tiling or penalties');
	}

	//recover elements distribution
	SIAflag  = FlagElements(md,options.getfieldvalue('SIA',''));
	SSAflag  = FlagElements(md,options.getfieldvalue('SSA',''));
	HOflag   = FlagElements(md,options.getfieldvalue('HO',''));
	L1L2flag = FlagElements(md,options.getfieldvalue('L1L2',''));
	MOLHOflag = FlagElements(md,options.getfieldvalue('MOLHO',''));
	FSflag   = FlagElements(md,options.getfieldvalue('FS',''));
	filltype = options.getfieldvalue('fill','none');
	options.displayunused();

	//Flag the elements that have not been flagged as filltype
	if (filltype === 'SIA'){
		for(var i=0;i<md.mesh.numberofelements;i++)if(!(SSAflag[i] | HOflag[i]))SIAflag[i]=1;
	}
	else if (filltype === 'SSA'){
		for(var i=0;i<md.mesh.numberofelements;i++)if(!(SIAflag[i] | HOflag[i] | FSflag[i]))SSAflag[i]=1;
	}
	else if (filltype === 'HO'){
		for(var i=0;i<md.mesh.numberofelements;i++)if(!(SIAflag[i] | SSAflag[i] | FSflag[i]))HOflag[i]=1;
	}

	//check that each element has at least one flag
	for(var i=0;i<md.mesh.numberofelements;i++)if((SIAflag[i] + SSAflag[i] + HOflag[i] + L1L2flag[i] + MOLHOflag[i] + FSflag[i])==0)
	throw Error("elements type not assigned, supported models are 'SIA','SSA','HO' and 'FS'");

	//check that each element has only one flag
	if (ArrayAnyAboveStrict(ArrayXPY(SIAflag,SSAflag,HOflag,L1L2flag,MOLHOflag),1)){
		console.log('setflowequation warning message: some elements have several types, higher order type is used for them')

		for(var i=0;i<md.mesh.numberofelements;i++){
			if(SIAflag[i] & SSAflag[i])SIAflag[i]=0;
			if(SIAflag[i] & HOflag[i])SIAflag[i]=0;
			if(SSAflag[i] & HOflag[i])SSAflag[i]=0;
		}
	}

	//check that L1L2 is not coupled to any other model for now
	if (ArrayAnyEqual(L1L2flag,1) & ArrayAnyEqual(ArrayOr(SIAflag,SSAflag,HOflag,FSflag),1)) throw Error('L1L2 cannot be coupled to any other model');
	if (ArrayAnyEqual(MOLHOflag,1) & ArrayAnyEqual(ArrayOr(SIAflag,SSAflag,HOflag,FSflag),1)) throw Error('MOLHO cannot be coupled to any other model');

	//Check that no HO or FS for 2d mesh
	if (md.mesh.domaintype() == '2Dhorizontal'){
		for(var i=0;i<FSflag.length;i++){
			if(FSflag[i] | HOflag[i]) throw Error('FS and HO elements not allowed in 2d mesh, extrude it first')
		}
	}

	//FS can only be used alone for now:
	if (ArrayAnyEqual(FSflag,1) & ArrayAnyEqual(SIAflag,1)) throw Error('FS cannot be used with any other model for now, put FS everywhere')

	//Initialize node fields
	nodeonSIA=NewArrayFill(md.mesh.numberofvertices,0);
	pos=ArrayFind(SIAflag,1);
	for(var i=0;i<pos.length;i++) for(var j=0;j<md.mesh.elements[0].length;j++) nodeonSIA[md.mesh.elements[pos[i]][j]-1]=1;
	
	nodeonSSA=NewArrayFill(md.mesh.numberofvertices,0);
	pos=ArrayFind(SSAflag,1);
	for(var i=0;i<pos.length;i++) for(var j=0;j<md.mesh.elements[0].length;j++) nodeonSSA[md.mesh.elements[pos[i]][j]-1]=1;
	
	nodeonHO=NewArrayFill(md.mesh.numberofvertices,0);
	pos=ArrayFind(HOflag,1);
	for(var i=0;i<pos.length;i++) for(var j=0;j<md.mesh.elements[0].length;j++) nodeonHO[md.mesh.elements[pos[i]][j]-1]=1;
	
	nodeonL1L2=NewArrayFill(md.mesh.numberofvertices,0);
	pos=ArrayFind(L1L2flag,1);
	for(var i=0;i<pos.length;i++) for(var j=0;j<md.mesh.elements[0].length;j++) nodeonL1L2[md.mesh.elements[pos[i]][j]-1]=1;

	nodeonMOLHO=NewArrayFill(md.mesh.numberofvertices,0);
	pos=ArrayFind(MOLHOflag,1);
	for(var i=0;i<pos.length;i++) for(var j=0;j<md.mesh.elements[0].length;j++) nodeonMOLHO[md.mesh.elements[pos[i]][j]-1]=1;

	nodeonFS=NewArrayFill(md.mesh.numberofvertices,0);
	noneflag=NewArrayFill(md.mesh.numberofvertices,0);
	
	
	//First modify FSflag to get rid of elements contrained everywhere (spc + border with HO or SSA)
	if (ArrayAnyEqual(FSflag,1)){
		throw Error("FS elements not supported yet!");
		/*fullspcnodes=double((~isnan(md.stressbalance.spcvx)+~isnan(md.stressbalance.spcvy)+~isnan(md.stressbalance.spcvz))==3 | (nodeonHO & nodeonFS));         //find all the nodes on the boundary of the domain without icefront
		fullspcelems=double(sum(fullspcnodes(md.mesh.elements),2)==6);         //find all the nodes on the boundary of the domain without icefront
		FSflag(find(fullspcelems))=0;
		nodeonFS(md.mesh.elements(find(FSflag),:))=1;*/
	}

	//Then complete with NoneApproximation or the other model used if there is no FS
	if (ArrayAnyEqual(FSflag,1)){
		throw Error("FS elements not supported yet!");
		/*if any(HOflag), //fill with HO
			HOflag(~FSflag)=1;
			nodeonHO(md.mesh.elements(find(HOflag),:))=1;
		elseif any(SSAflag), //fill with SSA
			SSAflag(~FSflag)=1;
			nodeonSSA(md.mesh.elements(find(SSAflag),:))=1;
		else //fill with none 
			noneflag(find(~FSflag))=1;
		end*/
	}

	//Now take care of the coupling between SSA and HO
	md.stressbalance.vertex_pairing=[];
	nodeonSSAHO=NewArrayFill(md.mesh.numberofvertices,0);
	nodeonHOFS=NewArrayFill(md.mesh.numberofvertices,0);
	nodeonSSAFS=NewArrayFill(md.mesh.numberofvertices,0);
	SSAHOflag=NewArrayFill(md.mesh.numberofelements,0);
	SSAFSflag=NewArrayFill(md.mesh.numberofelements,0);
	HOFSflag=NewArrayFill(md.mesh.numberofelements,0);

	/*if strcmpi(coupling_method,'penalties'),
		//Create the border nodes between HO and SSA and extrude them
		numnodes2d=md.mesh.numberofvertices2d;
		numlayers=md.mesh.numberoflayers;
		bordernodes2d=find(nodeonHO(1:numnodes2d) & nodeonSSA(1:numnodes2d)); //Nodes connected to two different types of elements

		//initialize and fill in penalties structure
		if ~isnan(bordernodes2d),
			penalties=[];
			for	i=1:numlayers-1,
				penalties=[penalties; [bordernodes2d bordernodes2d+md.mesh.numberofvertices2d*(i)]];
			end
			md.stressbalance.vertex_pairing=penalties;
		end
	elseif strcmpi(coupling_method,'tiling'),
		if any(SSAflag) & any(HOflag), //coupling SSA HO
			//Find node at the border
			nodeonSSAHO(find(nodeonSSA & nodeonHO))=1;
			//SSA elements in contact with this layer become SSAHO elements
			matrixelements=ismember(md.mesh.elements,find(nodeonSSAHO));
			commonelements=sum(matrixelements,2)~=0;
			commonelements(find(HOflag))=0; //only one layer: the elements previously in SSA
			SSAflag(find(commonelements))=0; //these elements are now SSAHOelements
			SSAHOflag(find(commonelements))=1;
			nodeonSSA(:)=0;
			nodeonSSA(md.mesh.elements(find(SSAflag),:))=1;

			//rule out elements that don't touch the 2 boundaries
			pos=find(SSAHOflag);
			elist=zeros(length(pos),1);
			elist = elist + any(sum(nodeonSSA(md.mesh.elements(pos,:)),2),2);
			elist = elist - any(sum(nodeonHO(md.mesh.elements(pos,:))  ,2),2);
			pos1=find(elist==1);
			SSAflag(pos(pos1))=1;
			SSAHOflag(pos(pos1))=0;
			pos2=find(elist==-1);
			HOflag(pos(pos2))=1;
			SSAHOflag(pos(pos2))=0;

			//Recompute nodes associated to these elements
			nodeonSSA(:)=0;
			nodeonSSA(md.mesh.elements(find(SSAflag),:))=1;
			nodeonHO(:)=0;
			nodeonHO(md.mesh.elements(find(HOflag),:))=1;
			nodeonSSAHO(:)=0;
			nodeonSSAHO(md.mesh.elements(find(SSAHOflag),:))=1;

		elseif any(HOflag) & any(FSflag), //coupling HO FS
			//Find node at the border
			nodeonHOFS(find(nodeonHO & nodeonFS))=1;
			//FS elements in contact with this layer become HOFS elements
			matrixelements=ismember(md.mesh.elements,find(nodeonHOFS));
			commonelements=sum(matrixelements,2)~=0;
			commonelements(find(HOflag))=0; //only one layer: the elements previously in SSA
			FSflag(find(commonelements))=0; //these elements are now SSAHOelements
			HOFSflag(find(commonelements))=1;
			nodeonFS=zeros(md.mesh.numberofvertices,1);
			nodeonFS(md.mesh.elements(find(FSflag),:))=1;

			//rule out elements that don't touch the 2 boundaries
			pos=find(HOFSflag);
			elist=zeros(length(pos),1);
			elist = elist + any(sum(nodeonFS(md.mesh.elements(pos,:)),2),2);
			elist = elist - any(sum(nodeonHO(md.mesh.elements(pos,:)),2),2);
			pos1=find(elist==1);
			FSflag(pos(pos1))=1;
			HOFSflag(pos(pos1))=0;
			pos2=find(elist==-1);
			HOflag(pos(pos2))=1;
			HOFSflag(pos(pos2))=0;

			//Recompute nodes associated to these elements
			nodeonFS(:)=0;
			nodeonFS(md.mesh.elements(find(FSflag),:))=1;
			nodeonHO(:)=0;
			nodeonHO(md.mesh.elements(find(HOflag),:))=1;
			nodeonHOFS(:)=0;
			nodeonHOFS(md.mesh.elements(find(HOFSflag),:))=1;

		elseif any(FSflag) & any(SSAflag),
			//Find node at the border
			nodeonSSAFS(find(nodeonSSA & nodeonFS))=1;
			//FS elements in contact with this layer become SSAFS elements
			matrixelements=ismember(md.mesh.elements,find(nodeonSSAFS));
			commonelements=sum(matrixelements,2)~=0;
			commonelements(find(SSAflag))=0; //only one layer: the elements previously in SSA
			FSflag(find(commonelements))=0; //these elements are now SSASSAelements
			SSAFSflag(find(commonelements))=1;
			nodeonFS=zeros(md.mesh.numberofvertices,1);
			nodeonFS(md.mesh.elements(find(FSflag),:))=1;

			//rule out elements that don't touch the 2 boundaries
			pos=find(SSAFSflag);
			elist=zeros(length(pos),1);
			elist = elist + any(sum(nodeonSSA(md.mesh.elements(pos,:)),2),2);
			elist = elist - any(sum(nodeonFS(md.mesh.elements(pos,:))  ,2),2);
			pos1=find(elist==1);
			SSAflag(pos(pos1))=1;
			SSAFSflag(pos(pos1))=0;
			pos2=find(elist==-1);
			FSflag(pos(pos2))=1;
			SSAFSflag(pos(pos2))=0;

			//Recompute nodes associated to these elements
			nodeonSSA(:)=0;
			nodeonSSA(md.mesh.elements(find(SSAflag),:))=1;
			nodeonFS(:)=0;
			nodeonFS(md.mesh.elements(find(FSflag),:))=1;
			nodeonSSAFS(:)=0;
			nodeonSSAFS(md.mesh.elements(find(SSAFSflag),:))=1;

		elseif any(FSflag) & any(SIAflag),
			error('type of coupling not supported yet');
		end
	end*/

	//Create element equations
	md.flowequation.element_equation=NewArrayFill(md.mesh.numberofelements,0);
	pos=ArrayFind(noneflag,1);for(var i=0;i<pos.length;i++)md.flowequation.element_equation[pos[i]]=0;
	pos=ArrayFind(SIAflag,1);for(var i=0;i<pos.length;i++)md.flowequation.element_equation[pos[i]]=1;
	pos=ArrayFind(SSAflag,1);for(var i=0;i<pos.length;i++)md.flowequation.element_equation[pos[i]]=2;
	pos=ArrayFind(L1L2flag,1);for(var i=0;i<pos.length;i++)md.flowequation.element_equation[pos[i]]=3;
	pos=ArrayFind(MOLHOflag,1);for(var i=0;i<pos.length;i++)md.flowequation.element_equation[pos[i]]=4;
	pos=ArrayFind(HOflag,1);for(var i=0;i<pos.length;i++)md.flowequation.element_equation[pos[i]]=5;
	pos=ArrayFind(FSflag,1);for(var i=0;i<pos.length;i++)md.flowequation.element_equation[pos[i]]=6;
	pos=ArrayFind(SSAHOflag,1);for(var i=0;i<pos.length;i++)md.flowequation.element_equation[pos[i]]=7;
	pos=ArrayFind(SSAFSflag,1);for(var i=0;i<pos.length;i++)md.flowequation.element_equation[pos[i]]=8;
	pos=ArrayFind(HOFSflag,1);for(var i=0;i<pos.length;i++)md.flowequation.element_equation[pos[i]]=9;


	//border
	md.flowequation.borderHO=nodeonHO;
	md.flowequation.borderSSA=nodeonSSA;
	md.flowequation.borderFS=nodeonFS;
	

	//Create vertices_type
	md.flowequation.vertex_equation=NewArrayFill(md.mesh.numberofvertices,0);

	pos=ArrayFind(nodeonSSA,1);for(var i=0;i<pos.length;i++)md.flowequation.vertex_equation[pos[i]]=2;
	pos=ArrayFind(nodeonL1L2,1);for(var i=0;i<pos.length;i++)md.flowequation.vertex_equation[pos[i]]=3;
	pos=ArrayFind(nodeonMOLHO,1);for(var i=0;i<pos.length;i++)md.flowequation.vertex_equation[pos[i]]=4;
	pos=ArrayFind(nodeonHO,1);for(var i=0;i<pos.length;i++)md.flowequation.vertex_equation[pos[i]]=5;
	pos=ArrayFind(nodeonFS,1);for(var i=0;i<pos.length;i++)md.flowequation.vertex_equation[pos[i]]=6;
	//DO SIA LAST! Otherwise spcs might not be set up correctly (SIA should have priority)
	pos=ArrayFind(nodeonSIA,1);for(var i=0;i<pos.length;i++)md.flowequation.vertex_equation[pos[i]]=1;
	if (ArrayAnyEqual(FSflag,1)){
		pos=ArrayFind(nodeonFS==0);
		if(ArrayAnyEqual(HOflag,0) & ArrayAnyEqual(SSA,0)){
			for(var i=0;i<pos.length;i++)md.flowequation.vertex_equation[pos[i]]=0;
		}
	}

	pos=ArrayFind(nodeonSSAHO,1);for(var i=0;i<pos.length;i++)md.flowequation.vertex_equation[pos[i]]=7;
	pos=ArrayFind(nodeonHOFS,1);for(var i=0;i<pos.length;i++)md.flowequation.vertex_equation[pos[i]]=8;
	pos=ArrayFind(nodeonSSAFS,2);for(var i=0;i<pos.length;i++)md.flowequation.vertex_equation[pos[i]]=9;

	//figure out solution types
	md.flowequation.isSIA  = ArrayAnyEqual(md.flowequation.element_equation,1);
	md.flowequation.isSSA  = ArrayAnyEqual(md.flowequation.element_equation,2);
	md.flowequation.isL1L2 = ArrayAnyEqual(md.flowequation.element_equation,3);
	md.flowequation.isMOLHO = ArrayAnyEqual(md.flowequation.element_equation,4);
	md.flowequation.isHO   = ArrayAnyEqual(md.flowequation.element_equation,5);
	md.flowequation.isFS   = ArrayAnyEqual(md.flowequation.element_equation,6);
	return

	//Check that tiling can work:
	/*if any(md.flowequation.borderSSA) & any(md.flowequation.borderHO) & any(md.flowequation.borderHO + md.flowequation.borderSSA ~=1),
		error('error coupling domain too irregular');
	end
	if any(md.flowequation.borderSSA) & any(md.flowequation.borderFS) & any(md.flowequation.borderFS + md.flowequation.borderSSA ~=1),
		error('error coupling domain too irregular');
	end
	if any(md.flowequation.borderFS) & any(md.flowequation.borderHO) & any(md.flowequation.borderHO + md.flowequation.borderFS~=1),
		error('error coupling domain too irregular');
	end*/
}
