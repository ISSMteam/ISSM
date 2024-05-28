function md=setflowequation(md,varargin)
%SETFLOWEQUATION - associate a solution type to each element
%
%   This routine works like plotmodel: it works with an even number of inputs
%   'SIA','SSA','L1L2','MOLHO','HO','FS' and 'fill' are the possible options
%   that must be followed by the corresponding exp file or flags list
%   It can either be a domain file (argus type, .exp extension), or an array of element flags. 
%   If user wants every element outside the domain to be 
%   setflowequationd, add '~' to the name of the domain file (ex: '~HO.exp');
%   an empty string '' will be considered as an empty domain
%   a string 'all' will be considered as the entire domain
%   You can specify the type of coupling, 'penalties' or 'tiling', to use with the input 'coupling'
%   NB: L1L2 and MOLHO cannot currently be coupled to any other ice flow model
%
%   Usage:
%      md=setflowequation(md,varargin)
%
%   Example:
%      md=setflowequation(md,'HO','HO.exp',fill','SIA','coupling','tiling');

%some checks on list of arguments
if ((nargin<2) | (nargout~=1)),
	error('setflowequation error message');
end

%Process options
options=pairoptions(varargin{:});
options=deleteduplicates(options,1);

%Find_out what kind of coupling to use
coupling_method=getfieldvalue(options,'coupling','tiling');
if (~strcmpi(coupling_method,'tiling') & ~strcmpi(coupling_method,'penalties')),
	error('coupling type can only be: tiling or penalties');
end

%recover elements distribution
SIAflag  = FlagElements(md,getfieldvalue(options,'SIA',''));
SSAflag  = FlagElements(md,getfieldvalue(options,'SSA',''));
HOflag   = FlagElements(md,getfieldvalue(options,'HO',''));
L1L2flag = FlagElements(md,getfieldvalue(options,'L1L2',''));
MOLHOflag = FlagElements(md,getfieldvalue(options,'MOLHO',''));
FSflag   = FlagElements(md,getfieldvalue(options,'FS',''));
filltype = getfieldvalue(options,'fill','none');
displayunused(options);

%Flag the elements that have not been flagged as filltype
if strcmpi(filltype,'SIA'),
	SIAflag(find(~(SSAflag | HOflag)))=1;
elseif strcmpi(filltype,'SSA'),
	SSAflag(find(~(SIAflag | HOflag | FSflag)))=1;
elseif strcmpi(filltype,'HO'),
	HOflag(find(~(SIAflag | SSAflag | FSflag)))=1;
end

%check that each element has at least one flag
if any(SIAflag+SSAflag+HOflag+L1L2flag+MOLHOflag+FSflag==0),
	error('elements type not assigned, supported models are ''SIA'',''SSA'',''HO'',''MOLHO'' and ''FS''')
end

%check that each element has only one flag
if any(SIAflag+SSAflag+HOflag+L1L2flag+MOLHOflag+FSflag>1),
	disp('setflowequation.m: Warning: some elements have several types, higher order type is used for them')
	SIAflag(find(SIAflag & SSAflag))=0;
	SIAflag(find(SIAflag & HOflag))=0;
	SSAflag(find(SSAflag & HOflag))=0;
end

%check that L1L2 is not coupled to any other model for now
if any(L1L2flag) & any(SIAflag | SSAflag | HOflag | FSflag)
	error('L1L2 cannot be coupled to any other model');
end
if any(MOLHOflag) & any(SIAflag | SSAflag | HOflag | FSflag)
	error('MOLHO cannot be coupled to any other model');
end

%Check that no HO or FS for 2d mesh
if strcmp(domaintype(md.mesh),'2Dhorizontal')
	if any(FSflag | HOflag)
		error('FS and HO elements not allowed in 2d mesh, extrude it first')
	end
end

%FS can only be used alone for now:
if any(FSflag) &any(SIAflag),
	error('FS cannot be used with any other model for now, put FS everywhere')
end

%Initialize node fields
nodeonSIA=zeros(md.mesh.numberofvertices,1);  nodeonSIA(md.mesh.elements(find(SIAflag),:))=1;
nodeonSSA=zeros(md.mesh.numberofvertices,1);  nodeonSSA(md.mesh.elements(find(SSAflag),:))=1;
nodeonHO=zeros(md.mesh.numberofvertices,1);   nodeonHO(md.mesh.elements(find(HOflag),:))=1;
nodeonL1L2=zeros(md.mesh.numberofvertices,1); nodeonL1L2(md.mesh.elements(find(L1L2flag),:))=1;
nodeonMOLHO=zeros(md.mesh.numberofvertices,1); nodeonMOLHO(md.mesh.elements(find(MOLHOflag),:))=1;
nodeonFS=zeros(md.mesh.numberofvertices,1);
noneflag=zeros(md.mesh.numberofelements,1);

%First modify FSflag to get rid of elements contrained everywhere (spc + border with HO or SSA)
if any(FSflag),
	fullspcnodes=double((~isnan(md.stressbalance.spcvx)+~isnan(md.stressbalance.spcvy)+~isnan(md.stressbalance.spcvz))==3 | (nodeonHO & nodeonFS));         %find all the nodes on the boundary of the domain without icefront
	fullspcelems=double(sum(fullspcnodes(md.mesh.elements),2)==6);         %find all the nodes on the boundary of the domain without icefront
	FSflag(find(fullspcelems))=0;
	nodeonFS(md.mesh.elements(find(FSflag),:))=1;
end

%Then complete with NoneApproximation or the other model used if there is no FS
if any(FSflag), 
	if any(HOflag), %fill with HO
		HOflag(~FSflag)=1;
		nodeonHO(md.mesh.elements(find(HOflag),:))=1;
	elseif any(SSAflag), %fill with SSA
		SSAflag(~FSflag)=1;
		nodeonSSA(md.mesh.elements(find(SSAflag),:))=1;
	else %fill with none 
		noneflag(find(~FSflag))=1;
	end
end

%Now take care of the coupling between SSA and HO
if strcmpi(coupling_method,'penalties'),
	md.stressbalance.vertex_pairing=[];
end
nodeonSSAHO=zeros(md.mesh.numberofvertices,1);
nodeonHOFS=zeros(md.mesh.numberofvertices,1);
nodeonSSAFS=zeros(md.mesh.numberofvertices,1);
SSAHOflag=zeros(md.mesh.numberofelements,1);
SSAFSflag=zeros(md.mesh.numberofelements,1);
HOFSflag=zeros(md.mesh.numberofelements,1);
if strcmpi(coupling_method,'penalties'),
	%Create the border nodes between HO and SSA and extrude them
	numnodes2d=md.mesh.numberofvertices2d;
	numlayers=md.mesh.numberoflayers;
	bordernodes2d=find(nodeonHO(1:numnodes2d) & nodeonSSA(1:numnodes2d)); %Nodes connected to two different types of elements

	%initialize and fill in penalties structure
	if ~isnan(bordernodes2d),
		penalties=[];
		for	i=1:numlayers-1,
			penalties=[penalties; [bordernodes2d bordernodes2d+md.mesh.numberofvertices2d*(i)]];
		end
		md.stressbalance.vertex_pairing=penalties;
	end
elseif strcmpi(coupling_method,'tiling'),
	if any(SSAflag) & any(HOflag), %coupling SSA HO
		%Find node at the border
		nodeonSSAHO(find(nodeonSSA & nodeonHO))=1;
		%SSA elements in contact with this layer become SSAHO elements
		matrixelements=ismember(md.mesh.elements,find(nodeonSSAHO));
		commonelements=sum(matrixelements,2)~=0;
		commonelements(find(HOflag))=0; %only one layer: the elements previously in SSA
		SSAflag(find(commonelements))=0; %these elements are now SSAHOelements
		SSAHOflag(find(commonelements))=1;
		nodeonSSA(:)=0;
		nodeonSSA(md.mesh.elements(find(SSAflag),:))=1;

		%rule out elements that don't touch the 2 boundaries
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

		%Recompute nodes associated to these elements
		nodeonSSA(:)=0;
		nodeonSSA(md.mesh.elements(find(SSAflag),:))=1;
		nodeonHO(:)=0;
		nodeonHO(md.mesh.elements(find(HOflag),:))=1;
		nodeonSSAHO(:)=0;
		nodeonSSAHO(md.mesh.elements(find(SSAHOflag),:))=1;

	elseif any(HOflag) & any(FSflag), %coupling HO FS
		%Find node at the border
		nodeonHOFS(find(nodeonHO & nodeonFS))=1;
		%FS elements in contact with this layer become HOFS elements
		matrixelements=ismember(md.mesh.elements,find(nodeonHOFS));
		commonelements=sum(matrixelements,2)~=0;
		commonelements(find(HOflag))=0; %only one layer: the elements previously in SSA
		FSflag(find(commonelements))=0; %these elements are now SSAHOelements
		HOFSflag(find(commonelements))=1;
		nodeonFS=zeros(md.mesh.numberofvertices,1);
		nodeonFS(md.mesh.elements(find(FSflag),:))=1;

		%rule out elements that don't touch the 2 boundaries
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

		%Recompute nodes associated to these elements
		nodeonFS(:)=0;
		nodeonFS(md.mesh.elements(find(FSflag),:))=1;
		nodeonHO(:)=0;
		nodeonHO(md.mesh.elements(find(HOflag),:))=1;
		nodeonHOFS(:)=0;
		nodeonHOFS(md.mesh.elements(find(HOFSflag),:))=1;

	elseif any(FSflag) & any(SSAflag),
		%Find node at the border
		nodeonSSAFS(find(nodeonSSA & nodeonFS))=1;
		%FS elements in contact with this layer become SSAFS elements
		matrixelements=ismember(md.mesh.elements,find(nodeonSSAFS));
		commonelements=sum(matrixelements,2)~=0;
		commonelements(find(SSAflag))=0; %only one layer: the elements previously in SSA
		FSflag(find(commonelements))=0; %these elements are now SSASSAelements
		SSAFSflag(find(commonelements))=1;
		nodeonFS=zeros(md.mesh.numberofvertices,1);
		nodeonFS(md.mesh.elements(find(FSflag),:))=1;

		%rule out elements that don't touch the 2 boundaries
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

		%Recompute nodes associated to these elements
		nodeonSSA(:)=0;
		nodeonSSA(md.mesh.elements(find(SSAflag),:))=1;
		nodeonFS(:)=0;
		nodeonFS(md.mesh.elements(find(FSflag),:))=1;
		nodeonSSAFS(:)=0;
		nodeonSSAFS(md.mesh.elements(find(SSAFSflag),:))=1;

	elseif any(FSflag) & any(SIAflag),
		error('type of coupling not supported yet');
	end
end

%Create element equations
md.flowequation.element_equation=zeros(md.mesh.numberofelements,1);
md.flowequation.element_equation(find(noneflag))=0;
md.flowequation.element_equation(find(SIAflag))=1;
md.flowequation.element_equation(find(SSAflag))=2;
md.flowequation.element_equation(find(L1L2flag))=3;
md.flowequation.element_equation(find(MOLHOflag))=4;
md.flowequation.element_equation(find(HOflag))=5;
md.flowequation.element_equation(find(FSflag))=6;
md.flowequation.element_equation(find(SSAHOflag))=7;
md.flowequation.element_equation(find(SSAFSflag))=8;
md.flowequation.element_equation(find(HOFSflag))=9;

%border
md.flowequation.borderHO=nodeonHO;
md.flowequation.borderSSA=nodeonSSA;
md.flowequation.borderFS=nodeonFS;

%Create vertices_type
md.flowequation.vertex_equation=zeros(md.mesh.numberofvertices,1);
pos=find(nodeonSSA);  md.flowequation.vertex_equation(pos)=2;
pos=find(nodeonL1L2); md.flowequation.vertex_equation(pos)=3;
pos=find(nodeonMOLHO); md.flowequation.vertex_equation(pos)=4;
pos=find(nodeonHO);   md.flowequation.vertex_equation(pos)=5;
pos=find(nodeonFS);   md.flowequation.vertex_equation(pos)=6;
%DO SIA LAST! Otherwise spcs might not be set up correctly (SIA should have priority)
pos=find(nodeonSIA);
md.flowequation.vertex_equation(pos)=1;
if any(FSflag),
	pos=find(~nodeonFS);
	if(~any(HOflag) & ~any(SSAflag)),
		md.flowequation.vertex_equation(pos)=0;
	end
end
pos=find(nodeonSSAHO);
md.flowequation.vertex_equation(pos)=7;
pos=find(nodeonHOFS);
md.flowequation.vertex_equation(pos)=8;
pos=find(nodeonSSAFS);
md.flowequation.vertex_equation(pos)=9;

%figure out solution types
md.flowequation.isSIA  = double(any(md.flowequation.element_equation == 1));
md.flowequation.isSSA  = double(any(md.flowequation.element_equation == 2));
md.flowequation.isL1L2 = double(any(md.flowequation.element_equation == 3));
md.flowequation.isMOLHO = double(any(md.flowequation.element_equation == 4));
md.flowequation.isHO   = double(any(md.flowequation.element_equation == 5));
md.flowequation.isFS   = double(any(md.flowequation.element_equation == 6));

return

%Check that tiling can work:
if any(md.flowequation.borderSSA) & any(md.flowequation.borderHO) & any(md.flowequation.borderHO + md.flowequation.borderSSA ~=1),
	error('error coupling domain too irregular');
end
if any(md.flowequation.borderSSA) & any(md.flowequation.borderFS) & any(md.flowequation.borderFS + md.flowequation.borderSSA ~=1),
	error('error coupling domain too irregular');
end
if any(md.flowequation.borderFS) & any(md.flowequation.borderHO) & any(md.flowequation.borderHO + md.flowequation.borderFS~=1),
	error('error coupling domain too irregular');
end
