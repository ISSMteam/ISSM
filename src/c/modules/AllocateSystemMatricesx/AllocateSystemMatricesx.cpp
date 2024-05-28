/*!\file AllocateSystemMatricesx
 * \brief retrieve vector from inputs in elements
 */

#include "./AllocateSystemMatricesx.h"
#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"

void AllocateSystemMatricesx(Matrix<IssmDouble>** pKff,Matrix<IssmDouble>** pKfs,Vector<IssmDouble>** pdf,Vector<IssmDouble>** ppf,FemModel* femmodel){

	/*Intermediary*/
	int  fsize,ssize,flocalsize,slocalsize;
	int  connectivity, numberofdofspernode;
	int  m,n,M,N;
	int *d_nnz = NULL;
	int *o_nnz = NULL;

	/*output*/
	Matrix<IssmDouble> *Kff  = NULL;
	Matrix<IssmDouble> *Kfs  = NULL;
	Vector<IssmDouble> *pf   = NULL;
	Vector<IssmDouble> *df   = NULL;

	bool oldalloc=false;
	char* toolkittype=NULL;

	/*retrieve parameters: */
	femmodel->parameters->FindParam(&connectivity,MeshAverageVertexConnectivityEnum);

	/*retrieve node info*/
	fsize      = femmodel->nodes->NumberOfDofs(FsetEnum);
	ssize      = femmodel->nodes->NumberOfDofs(SsetEnum);
	flocalsize = femmodel->nodes->NumberOfDofsLocal(FsetEnum);
	slocalsize = femmodel->nodes->NumberOfDofsLocal(SsetEnum);

	numberofdofspernode=femmodel->nodes->MaxNumDofs(GsetEnum);

	/*if our matrices are coming from issm, we don't do dynamic allocation like Petsc
	 * does, and this routine is essentially useless. Force standard alloc in this case: */
	toolkittype=ToolkitOptions::GetToolkitType();

	if(oldalloc){
		if(pKff) Kff=new Matrix<IssmDouble>(fsize,fsize,connectivity,numberofdofspernode);
		if(pKfs) Kfs=new Matrix<IssmDouble>(fsize,ssize,connectivity,numberofdofspernode);
		if(pdf)  df =new Vector<IssmDouble>(fsize);
		if(ppf)  pf =new Vector<IssmDouble>(fsize);
	}
	else{
		if(pKff){
			m=flocalsize; n=flocalsize; /*local  sizes*/
			M=fsize;      N=fsize;      /*global sizes*/
			if(strcmp(toolkittype,"issm")==0){
				Kff=new Matrix<IssmDouble>(m,n,M,N,NULL,NULL);
			}
			else{
				MatrixNonzeros(&d_nnz,&o_nnz,femmodel,FsetEnum,FsetEnum);
				Kff=new Matrix<IssmDouble>(m,n,M,N,d_nnz,o_nnz);
				xDelete<int>(d_nnz);
				xDelete<int>(o_nnz);
			}
		}
		if(pKfs){
			m=flocalsize; n=slocalsize; /*local  sizes*/
			M=fsize;      N=ssize;      /*global sizes*/
			if(strcmp(toolkittype,"issm")==0){
				Kfs=new Matrix<IssmDouble>(m,n,M,N,NULL,NULL);
			}
			else{
				MatrixNonzeros(&d_nnz,&o_nnz,femmodel,FsetEnum,SsetEnum);
				Kfs=new Matrix<IssmDouble>(m,n,M,N,d_nnz,o_nnz);
				xDelete<int>(d_nnz);
				xDelete<int>(o_nnz);
			}
		}
		if(pdf) df =new Vector<IssmDouble>(flocalsize,fsize);
		if(ppf) pf =new Vector<IssmDouble>(flocalsize,fsize);
	}

	/*Free resources: */
	xDelete<char>(toolkittype);

	/*Allocate output pointers*/
	if(pKff) *pKff = Kff;
	if(pKfs) *pKfs = Kfs;
	if(pdf)  *pdf  = df;
	if(ppf)  *ppf  = pf;
}

void MatrixNonzeros(int** pd_nnz,int** po_nnz,FemModel* femmodel,int set1enum,int set2enum){

	/*Intermediary*/
	int      i,j,k,index,offset,count;
	int      d_nz,o_nz;
	Element *element            = NULL;
	Load    *load               = NULL;
	int     *head_e             = NULL;
	int     *next_e             = NULL;
	int     *count2offset_e     = NULL;
	int     *head_l             = NULL;
	int     *next_l             = NULL;
	int     *count2offset_l     = NULL;
	int     *lidlist            = NULL;

	/*output*/
	int *d_nnz = NULL;
	int *o_nnz = NULL;

	/*Get vector size and number of nodes*/
	int numnodes            = femmodel->nodes->NumberOfNodes();
	int localmasters        = femmodel->nodes->NumberOfNodesLocal();
	int localnumnodes       = femmodel->nodes->Size();
	int numberofdofspernode = femmodel->nodes->MaxNumDofs(GsetEnum);
	//int M                   = femmodel->nodes->NumberOfDofs(set1enum);
	int N                   = femmodel->nodes->NumberOfDofs(set2enum);
	int m                   = femmodel->nodes->NumberOfDofsLocal(set1enum);
	int n                   = femmodel->nodes->NumberOfDofsLocal(set2enum);
	int numnodesperelement  = femmodel->elements->MaxNumNodes();
	int numnodesperload     = femmodel->loads->MaxNumNodes();
	int elementssize        = femmodel->elements->Size();
	int loadssize           = femmodel->loads->Size();

	/*First, we are building chaining vectors so that we know what nodes are
	 * connected to what elements. These vectors are such that:
	 *   for(int i=head[id];i!=-1;i=next[i])
	 * will loop over all the elements that are connected to the node number
	 * id*/
	head_e         = xNew<int>(localnumnodes); for(i=0;i<localnumnodes;i++) head_e[i]=-1;
	next_e         = xNew<int>(elementssize*numnodesperelement);
	count2offset_e = xNew<int>(elementssize*numnodesperelement);
	k=0;
	i=-1;
	for(Object* & object : femmodel->elements->objects){
      element = xDynamicCast<Element*>(object);
		i+=1;
		int elementnumnodes = element->GetNumberOfNodes();
		lidlist = xNew<int>(elementnumnodes);
		element->GetNodesLidList(lidlist);

		for(j=0;j<elementnumnodes;j++){
			index = lidlist[j];
			_assert_(index>=0 && index<numnodes);

			count2offset_e[k]=i;
			next_e[k]=head_e[index];
			head_e[index]=k++;
		}
		k = k + (numnodesperelement-elementnumnodes);

		xDelete<int>(lidlist);
	}
	/*Chain for loads*/
	head_l         = xNew<int>(localnumnodes); for(i=0;i<localnumnodes;i++) head_l[i]=-1;
	next_l         = xNew<int>(loadssize*numnodesperload);
	count2offset_l = xNew<int>(loadssize*numnodesperload);
	k=0;
	i=-1;
	for(Object* & object : femmodel->loads->objects){
      load = xDynamicCast<Load*>(object);
		i+=1;
		int loadnumnodes = load->GetNumberOfNodes();
		lidlist = xNew<int>(loadnumnodes);
		load->GetNodesLidList(lidlist);

		for(j=0;j<loadnumnodes;j++){
			index = lidlist[j];
			_assert_(index>=0 && index<numnodes);

			count2offset_l[k]=i;
			next_l[k]=head_l[index];
			head_l[index]=k++;
		}
		k = k + (numnodesperload-loadnumnodes);

		xDelete<int>(lidlist);
	}

	/*OK now count number of dofs and flag each nodes for each node i*/
	bool *flags                = xNew<bool>(localnumnodes);
	int  *flagsindices         = xNew<int>(localnumnodes);
	int  *d_connectivity       = xNewZeroInit<int>(localnumnodes);
	int  *o_connectivity       = xNewZeroInit<int>(localnumnodes);
	int   flagsindices_counter;
	int   analysis_type;

	Vector<IssmDouble> *connectivity_clone= new Vector<IssmDouble>(localmasters,numnodes);
	femmodel->parameters->FindParam(&analysis_type,AnalysisTypeEnum);
	/*Resetting flags to false at each iteration takes a lot of time, so we keep track of the flags
	 * to reset in flagsindices, initialized with -1*/
	for(i = 0;i<localnumnodes;i++) flags[i]        = false;
	for(i = 0;i<localnumnodes;i++) flagsindices[i] = -1;

	/*Create connectivity vector*/
	for(Object* & object : femmodel->nodes->objects){
      Node* node = xDynamicCast<Node*>(object);
		int   lid = node->Lid();
		int   pid = node->Pid();
		/*Reinitialize flags to false*/
		j=0;
		while(j<localnumnodes){
			if(flagsindices[j]>=0){
				flags[flagsindices[j]] = false;
				flagsindices[j]        = -1;
				j++;
			}
			else{
				break;
			}
		}
		flagsindices_counter = 0;
		for(j=head_e[lid];j!=-1;j=next_e[j]){
			offset=count2offset_e[j];
			element=xDynamicCast<Element*>(femmodel->elements->GetObjectByOffset(offset));
			element->SetwiseNodeConnectivity(&d_nz,&o_nz,node,flags,flagsindices,&flagsindices_counter,set1enum,set2enum,analysis_type);
			if(node->IsClone()){
				connectivity_clone->SetValue(pid,d_nz+o_nz,ADD_VAL);
			}
			else{
				d_connectivity[lid]+=d_nz;
				o_connectivity[lid]+=o_nz;
			}
		}
		for(j=head_l[lid];j!=-1;j=next_l[j]){
			offset=count2offset_l[j];
			load=xDynamicCast<Load*>(femmodel->loads->GetObjectByOffset(offset));
			load->SetwiseNodeConnectivity(&d_nz,&o_nz,node,flags,flagsindices,&flagsindices_counter,set1enum,set2enum);
			if(node->IsClone()){
				connectivity_clone->SetValue(pid,d_nz+o_nz,ADD_VAL);
			}
			else{
				d_connectivity[lid]+=d_nz;
				o_connectivity[lid]+=o_nz;
			}
		}
	}
	xDelete<bool>(flags);
	xDelete<int>(flagsindices);
	xDelete<int>(count2offset_e);
	xDelete<int>(head_e);
	xDelete<int>(next_e);
	xDelete<int>(count2offset_l);
	xDelete<int>(head_l);
	xDelete<int>(next_l);

	/*sum over all cpus*/
	connectivity_clone->Assemble();
	IssmDouble* serial_connectivity_clone=NULL;
	femmodel->GetLocalVectorWithClonesNodes(&serial_connectivity_clone,connectivity_clone);
	delete connectivity_clone;

	if(set1enum==FsetEnum){
		count=0;
		d_nnz=xNew<int>(m);
		o_nnz=xNew<int>(m);
		for(Object* & object : femmodel->nodes->objects){
			Node* node = xDynamicCast<Node*>(object);
			int   lid = node->Lid();
			if(!node->IsClone()){
				int node_fsize = node->FSize();
				for(j=0;j<node_fsize;j++){
					_assert_(count<m);
					d_nnz[count]=numberofdofspernode*(d_connectivity[lid] + reCast<int>(serial_connectivity_clone[lid]));
					o_nnz[count]=numberofdofspernode*(o_connectivity[lid] + reCast<int>(serial_connectivity_clone[lid]));
					if(d_nnz[count]>n)   d_nnz[count]=n;
					if(o_nnz[count]>N-n) o_nnz[count]=N-n;
					count++;
				}
			}
		}
		_assert_(m==count);
	}
	else{
		_error_("STOP not implemented");
	}
	xDelete<int>(d_connectivity);
	xDelete<int>(o_connectivity);
	xDelete<IssmDouble>(serial_connectivity_clone);

	/*Allocate ouptput pointer*/
	*pd_nnz=d_nnz;
	*po_nnz=o_nnz;
}
