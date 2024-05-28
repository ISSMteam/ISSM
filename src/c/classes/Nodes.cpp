/*
 * \file Nodes.cpp
 * \brief: Implementation of Nodes class, derived from DataSet class.
 */

/*Headers*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif
#include "../shared/io/Comm/IssmComm.h"
#include "./Nodes.h"
#include "./Node.h"
using namespace std;

/*Object constructors and destructor*/
Nodes::Nodes(){/*{{{*/
	this->enum_type             = NodesEnum;
	this->common_recv           = NULL;
	this->common_recv_ids       = NULL;
	this->common_send           = NULL;
	this->common_send_ids       = NULL;
	this->numberofnodes         = -1;
	this->numberofnodes_local   = -1;
	this->numberofmasters_local = -1;
	return;
}
/*}}}*/
Nodes::~Nodes(){/*{{{*/
	int num_proc=IssmComm::GetSize();

	if(this->common_recv) xDelete<int>(common_recv);
	if(this->common_send) xDelete<int>(common_send);
	if(this->common_recv_ids){
		for(int i=0;i<num_proc;i++) if(common_recv_ids[i]) xDelete<int>(common_recv_ids[i]);
		xDelete<int*>(common_recv_ids);
	}
	if(this->common_send_ids){
		for(int i=0;i<num_proc;i++) if(common_send_ids[i]) xDelete<int>(common_send_ids[i]);
		xDelete<int*>(common_send_ids);
	}
	return;
}
/*}}}*/

/*Numerics*/
Nodes* Nodes::Copy() {/*{{{*/

	int num_proc = IssmComm::GetSize();

	/*Copy dataset*/
	Nodes* output=new Nodes();
	output->sorted=this->sorted;
	output->numsorted=this->numsorted;
	output->presorted=this->presorted;
	for(vector<Object*>::iterator obj=this->objects.begin() ; obj < this->objects.end(); obj++ ){
		output->AddObject((*obj)->copy());
	}

	/*Build id_offsets and sorted_ids*/
	int objsize = this->numsorted;
	output->id_offsets=NULL;
	output->sorted_ids=NULL;
	if(this->sorted && objsize>0 && this->id_offsets){
		output->id_offsets=xNew<int>(objsize);
		xMemCpy<int>(output->id_offsets,this->id_offsets,objsize);
	}
	if(this->sorted && objsize>0 && this->sorted_ids){
		output->sorted_ids=xNew<int>(objsize);
		xMemCpy<int>(output->sorted_ids,this->sorted_ids,objsize);
	}

	/*Copy other fields*/
	output->numberofnodes         = this->numberofnodes;
	output->numberofnodes_local   = this->numberofnodes_local;
	output->numberofmasters_local = this->numberofmasters_local;

	if(this->common_recv){
		output->common_recv=xNew<int>(num_proc);
		for(int i=0;i<num_proc;i++) output->common_recv[i]=this->common_recv[i];
	}
	if(this->common_send){
		output->common_send=xNew<int>(num_proc);
		for(int i=0;i<num_proc;i++) output->common_send[i]=this->common_send[i];
	}
	if(this->common_recv_ids){
		output->common_recv_ids = xNew<int*>(num_proc);
		for(int i=0;i<num_proc;i++){
			output->common_recv_ids[i]=xNew<int>(this->common_recv[i]);
			for(int j=0;j<this->common_recv[i];j++) output->common_recv_ids[i][j]=this->common_recv_ids[i][j];
		}
	}
	if(this->common_send_ids){
		output->common_send_ids = xNew<int*>(num_proc);
		for(int i=0;i<num_proc;i++){
			output->common_send_ids[i]=xNew<int>(this->common_send[i]);
			for(int j=0;j<this->common_send[i];j++) output->common_send_ids[i][j]=this->common_send_ids[i][j];
		}
	}

	return output;
}
/*}}}*/
void  Nodes::Marshall(MarshallHandle* marshallhandle){ /*{{{*/

	int object_enum = NodesEnum;
	marshallhandle->call(object_enum);

	marshallhandle->call(numberofnodes);
	marshallhandle->call(numberofnodes_local);
	marshallhandle->call(numberofmasters_local);

	/*Check that restart is compatible!*/
	int num_procs=IssmComm::GetSize();
	int test = num_procs;
	marshallhandle->call(test);
	if(test!=num_procs) _error_("number of cores is not the same as before");

	DataSet::Marshall(marshallhandle);

	if(marshallhandle->OperationNumber() == MARSHALLING_LOAD){
		this->common_recv_ids = xNew<int*>(num_procs);
		this->common_send_ids = xNew<int*>(num_procs);
		for(int i=0;i<num_procs;i++){
			this->common_recv_ids[i] = NULL;
			this->common_send_ids[i] = NULL;
		}
	}

	/*Stop here if no nodes*/
	if(this->Size()==0) return;

	marshallhandle->call(this->common_recv,num_procs);
	marshallhandle->call(this->common_send,num_procs);
	for(int i=0;i<num_procs;i++){
		if(this->common_recv[i]) marshallhandle->call(this->common_recv_ids[i],this->common_recv[i]);
		if(this->common_send[i]) marshallhandle->call(this->common_send_ids[i],this->common_send[i]);
	}
}
/*}}}*/
void  Nodes::DistributeDofs(int setenum){/*{{{*/

	/*recover my_rank:*/
	ISSM_MPI_Status status;
	int my_rank   = IssmComm::GetRank();
	int num_procs = IssmComm::GetSize();

	/*Now, Build local dofs for masters first*/
	int  dofcount=0;
	for(Object* & object : this->objects){
      Node* node = xDynamicCast<Node*>(object);
		if(!node->IsClone()) node->DistributeLocalDofs(&dofcount,setenum);
	}
	/*Build local dofs for clones, they always will be at the end*/
	int dofcount_local = dofcount;
	for(Object* & object : this->objects){
		Node* node = xDynamicCast<Node*>(object);
		if(node->IsClone()) node->DistributeLocalDofs(&dofcount_local,setenum);
	}

	/* Now every object has distributed dofs, but locally, and with a dof count starting from
	 * 0. This means the dofs between all the cpus are not unique. We now offset the dofs of each
	 * cpus by the total last (master) dofs of the previus cpu, starting from 0.
	 * First: get number of dofs for each cpu*/
	int* alldofcount=xNew<int>(num_procs);
	ISSM_MPI_Gather(&dofcount,1,ISSM_MPI_INT,alldofcount,1,ISSM_MPI_INT,0,IssmComm::GetComm());
	ISSM_MPI_Bcast(alldofcount,num_procs,ISSM_MPI_INT,0,IssmComm::GetComm());

	/* Every cpu should start its own dof count at the end of the dofcount from cpu-1*/
	int offset=0;
	for(int i=0;i<my_rank;i++) offset+=alldofcount[i];
	xDelete<int>(alldofcount);

	for(Object* & object : this->objects){
		Node* node = xDynamicCast<Node*>(object);
		node->DistributeGlobalDofsMasters(offset,setenum);
	}

	/* Finally, remember that cpus may have skipped some objects, because they were clones. For every
	 * object that is not a clone, tell them to show their dofs, so that later on, they can get picked
	 * up by their clones: */
	int maxdofspernode  = this->MaxNumDofs(GsetEnum);
	int **send_truedofs = xNewZeroInit<int*>(num_procs);
	int  *recv_truedofs = xNewZeroInit<int>(this->Size()*maxdofspernode);
	ISSM_MPI_Request  *send_requests = xNew<ISSM_MPI_Request>(num_procs);
	for (int rank = 0;rank<num_procs;rank++) send_requests[rank] = ISSM_MPI_REQUEST_NULL;

	for(int rank=0;rank<num_procs;rank++){
		if(this->common_send[rank]){
			int  numids = this->common_send[rank];
			send_truedofs[rank] = xNew<int>(numids*maxdofspernode);
			for(int i=0;i<numids;i++){
				Node* node=xDynamicCast<Node*>(this->GetObjectByOffset(this->common_send_ids[rank][i]));
				node->ShowMasterDofs(&send_truedofs[rank][i*maxdofspernode+0],setenum);
			}
			ISSM_MPI_Isend(send_truedofs[rank],numids*maxdofspernode,ISSM_MPI_INT,rank,0,IssmComm::GetComm(),&send_requests[rank]);
		}
	}
	for(int rank=0;rank<num_procs;rank++){
		if(this->common_recv[rank]){
			int  numids = this->common_recv[rank];
			ISSM_MPI_Recv(recv_truedofs,numids*maxdofspernode,ISSM_MPI_INT,rank,0,IssmComm::GetComm(),&status);
			for(int i=0;i<numids;i++){
				Node* node=xDynamicCast<Node*>(this->GetObjectByOffset(this->common_recv_ids[rank][i]));
				node->UpdateCloneDofs(&recv_truedofs[i*maxdofspernode+0],setenum);
			}
		}
	}
	xDelete<int>(recv_truedofs);
	for(int rank=0;rank<num_procs;rank++){
		if(this->common_send[rank]) ISSM_MPI_Wait(&send_requests[rank],&status);
		xDelete<int>(send_truedofs[rank]);
	}
	xDelete<int*>(send_truedofs);
	xDelete<ISSM_MPI_Request>(send_requests);

	/*Update indexingupdateflag*/
	for(Object* & object : this->objects){
		Node* node = xDynamicCast<Node*>(object);
		node->ReindexingDone();
	}
}
/*}}}*/
void  Nodes::Finalize(){/*{{{*/

	/*Here we do 4 things:
	 * - count all nodes once for all so that we do not need to call MPI
	 *   every time we need to know the total number of vertices
	 * - Distribute lids (local ids): masters first, slaves second
	 * - Distribute pids (parallel ids)
	 * - Distribute Gset once for all
	 */

	/*recover my_rank:*/
	ISSM_MPI_Status status;
	int my_rank   = IssmComm::GetRank();
	int num_procs = IssmComm::GetSize();

	/*1. set number of nodes once for all*/
	this->numberofnodes_local=this->Size();
	this->numberofmasters_local=0;
	for(Object* & object : this->objects){
		Node* node = xDynamicCast<Node*>(object);
		if(!node->clone) this->numberofmasters_local++;
	}
	ISSM_MPI_Allreduce((void*)&this->numberofmasters_local,(void*)&this->numberofnodes,1,ISSM_MPI_INT,ISSM_MPI_SUM,IssmComm::GetComm());

	/*2. Distribute lids (First: masters, then clones)*/
	int lid = 0;
	for(Object* & object : this->objects){
		Node* node = xDynamicCast<Node*>(object);
		if(!node->clone) node->lid=lid++;
	}
	for(Object* & object : this->objects){
		Node* node = xDynamicCast<Node*>(object);
		if(node->clone) node->lid=lid++;
	}

	/*3. Distribute pids based on lids and offsets*/
	int* all_num_masters=xNew<int>(num_procs);
	ISSM_MPI_Gather(&this->numberofmasters_local,1,ISSM_MPI_INT,all_num_masters,1,ISSM_MPI_INT,0,IssmComm::GetComm());
	ISSM_MPI_Bcast(all_num_masters,num_procs,ISSM_MPI_INT,0,IssmComm::GetComm());
	int offset=0;
	for(int i=0;i<my_rank;i++) offset+=all_num_masters[i];
	xDelete<int>(all_num_masters);

	for(Object* & object : this->objects){
		Node* node = xDynamicCast<Node*>(object);
		node->pid = node->lid+offset;
	}

	/* Share pids of masters and update pids of clones*/
	int **send_truepids = xNewZeroInit<int*>(num_procs);
	int  *recv_truepids = xNewZeroInit<int>(this->Size());
	ISSM_MPI_Request* send_requests = xNew<ISSM_MPI_Request>(num_procs);
	for(int rank=0;rank<num_procs;rank++) send_requests[rank] = ISSM_MPI_REQUEST_NULL;
	for(int rank=0;rank<num_procs;rank++){
		if(this->common_send[rank]){
			int  numids = this->common_send[rank];
			send_truepids[rank] = xNew<int>(numids);
			for(int i=0;i<numids;i++){
				Node* node=xDynamicCast<Node*>(this->GetObjectByOffset(this->common_send_ids[rank][i]));
				send_truepids[rank][i] = node->pid;
			}
			ISSM_MPI_Isend(send_truepids[rank],numids,ISSM_MPI_INT,rank,0,IssmComm::GetComm(),&send_requests[rank]);
		}
	}
	for(int rank=0;rank<num_procs;rank++){
		if(this->common_recv[rank]){
			int  numids = this->common_recv[rank];
			ISSM_MPI_Recv(recv_truepids,numids,ISSM_MPI_INT,rank,0,IssmComm::GetComm(),&status);
			for(int i=0;i<numids;i++){
				Node* node=xDynamicCast<Node*>(this->GetObjectByOffset(this->common_recv_ids[rank][i]));
				node->pid = recv_truepids[i];
			}
		}
	}
	xDelete<int>(recv_truepids);
	for(int rank=0;rank<num_procs;rank++){
		if(this->common_send[rank]) ISSM_MPI_Wait(&send_requests[rank],&status);
		xDelete<int>(send_truepids[rank]);
	}
	xDelete<int*>(send_truepids);
	xDelete<ISSM_MPI_Request>(send_requests);

	/*4. Distribute G dofs once for all*/
	//this->DistributeDofs(GsetEnum);

	return;
}/*}}}*/
int   Nodes::MaxNumDofs(int setenum){/*{{{*/

	int max=0;
	int allmax;

	/*Now go through all nodes, and get how many dofs they own, unless they are clone nodes: */
	for(Object* & object : this->objects){
		Node* node = xDynamicCast<Node*>(object);

		int numdofs=node->GetNumberOfDofs(NoneApproximationEnum,setenum);
		if(numdofs>max) max=numdofs;
	}

	/*Grab max of all cpus: */
	ISSM_MPI_Allreduce((void*)&max,(void*)&allmax,1,ISSM_MPI_INT,ISSM_MPI_MAX,IssmComm::GetComm());
	max=allmax;

	return max;
}
/*}}}*/
int   Nodes::NumberOfDofs(int setenum){/*{{{*/

	int   allnumdofs;

	/*Get number of dofs on current cpu (excluding clones)*/
	int numdofs=this->NumberOfDofsLocal(setenum);

	/*Gather from all cpus: */
	ISSM_MPI_Allreduce ((void*)&numdofs,(void*)&allnumdofs,1,ISSM_MPI_INT,ISSM_MPI_SUM,IssmComm::GetComm());
	return allnumdofs;
}
/*}}}*/
int   Nodes::NumberOfDofsLocal(int setenum){/*{{{*/

	int   numdofs=0;

	/*Now go through all nodes, and get how many dofs they own, unless they are clone nodes: */
	for(Object* & object : this->objects){
		Node* node = xDynamicCast<Node*>(object);

		/*Ok, this object is a node, ask it to plug values into partition: */
		if (!node->IsClone()){
			numdofs+=node->GetNumberOfDofs(NoneApproximationEnum,setenum);
		}
	}

	return numdofs;
}
/*}}}*/
int   Nodes::NumberOfDofsLocalAll(int setenum){/*{{{*/

	/*go through all nodes, and get how many dofs they own, unless they are clone nodes: */
	int numdofs=0;
	for(Object* & object : this->objects){
		Node* node = xDynamicCast<Node*>(object);
		numdofs+=node->GetNumberOfDofs(NoneApproximationEnum,setenum);
	}
	return numdofs;
}
/*}}}*/
int   Nodes::NumberOfNodes(void){/*{{{*/

	return this->numberofnodes;
}
/*}}}*/
int   Nodes::NumberOfNodesLocal(void){/*{{{*/

	return this->numberofmasters_local;
}
/*}}}*/
int   Nodes::NumberOfNodesLocalAll(void){/*{{{*/

	return this->numberofnodes_local;
}
/*}}}*/
bool  Nodes::RequiresDofReindexing(void){/*{{{*/

	int flag = 0;
	int allflag;

	/*Now go through all nodes, and get how many dofs they own, unless they are clone nodes: */
	for(Object* & object : this->objects){
		Node* node = xDynamicCast<Node*>(object);
		if(node->RequiresDofReindexing()){
			flag = 1;
			break;
		}
	}

	/*Grab max of all cpus: */
	ISSM_MPI_Allreduce((void*)&flag,(void*)&allflag,1,ISSM_MPI_INT,ISSM_MPI_MAX,IssmComm::GetComm());

	if(allflag) return true;
	else        return false;
}
/*}}}*/

void  Nodes::CheckDofListAcrossPartitions(void){/*{{{*/

	/*recover my_rank:*/
	ISSM_MPI_Status status;
	int my_rank   = IssmComm::GetRank();
	int num_procs = IssmComm::GetSize();

	/*Display message*/
	if(VerboseModule()) _printf0_("   Checking degrees of freedom across partitions\n");

	/*Allocate vector to check degrees of freedom*/
	int gsize      = this->NumberOfDofs(GsetEnum);
	int glocalsize = this->NumberOfDofsLocal(GsetEnum);
	Vector<IssmDouble>* dofs_check=new Vector<IssmDouble>(glocalsize,gsize);

	/*First, go over all nodes, and masters can write their f dof and -1 for s-set*/
	for(Object* & object : this->objects){
		Node* node = xDynamicCast<Node*>(object);

		/*Skip if clone (will check later)*/
		if(node->IsClone()) continue;

		/*Write degree of freedom if active*/
		int count = 0;
		for(int j=0;j<node->gsize;j++){
			if(node->f_set[j]){
				if(node->s_set[j]) _error_("a degree of freedom is both in f and s set!");
				dofs_check->SetValue(node->gdoflist[j],reCast<IssmDouble>(node->fdoflist[count]),INS_VAL);
				count++;
			}
			else{
				if(node->s_set[j]==0) _error_("a degree of freedom is neither in f nor in s set!");
				dofs_check->SetValue(node->gdoflist[j],-1.,INS_VAL);
			}
		}
	}
	dofs_check->Assemble();

	/*Get local vector with both masters and slaves:*/
	IssmDouble *local_dofs_check = NULL;
	this->GetLocalVectorWithClonesGset(&local_dofs_check,dofs_check);
	delete dofs_check;

	/*Second, go over all nodes, and check that we still have what's expected...*/
	for(Object* & object : this->objects){
		Node* node = xDynamicCast<Node*>(object);

		/*Write degree of freedom if active*/
		int countg = 0;
		int countf = 0;
		int counts = 0;
		for(int j=0;j<node->gsize;j++){
			int index = node->gdoflist_local[countg];
			if(node->f_set[j]){
				if(reCast<int>(local_dofs_check[index]) != node->fdoflist[countf]){
					_error_("Dof #"<<j<<" of node sid "<<node->Sid()<<" not consistent: "<<local_dofs_check[index]<<"!="<<node->fdoflist[countf]);
				}
				countf++;
			}
			else{
				if(local_dofs_check[index] != -1.){
					_error_("Dof #"<<j<<" of node sid "<<node->Sid()<<" not consistently in s set");
				}
				counts++;
			}
			countg++;
		}
	}

	/*cleanup and return*/
	xDelete<IssmDouble>(local_dofs_check);
}/*}}}*/
void  Nodes::GetLocalVectorWithClonesGset(IssmDouble** plocal_ug,Vector<IssmDouble> *ug){/*{{{*/

	/*recover my_rank:*/
	ISSM_MPI_Status status;
	int my_rank   = IssmComm::GetRank();
	int num_procs = IssmComm::GetSize();

	/*retrieve node info*/
	int glocalsize         = this->NumberOfDofsLocalAll(GsetEnum);
	int glocalsize_masters = this->NumberOfDofsLocal(GsetEnum);
	int maxdofspernode     = this->MaxNumDofs(GsetEnum);

	/*Get local vector of ug*/
	int        *indices_ug_masters = NULL;
	IssmDouble *local_ug_masters   = NULL;
	ug->GetLocalVector(&local_ug_masters,&indices_ug_masters);
	_assert_(glocalsize_masters==indices_ug_masters[glocalsize_masters-1] - indices_ug_masters[0]+1);
	xDelete<int>(indices_ug_masters);

	/*Now, extend vectors to account for clones (make vectors longer, for clones at the end)*/
	IssmDouble *local_ug  = xNew<IssmDouble>(glocalsize);
	xMemCpy<IssmDouble>(local_ug,local_ug_masters,glocalsize_masters);
	xDelete<IssmDouble>(local_ug_masters);

	/*Now send and receive ug for nodes on partition edge*/
	IssmDouble **send_buffers = xNewZeroInit<IssmDouble*>(num_procs);
	IssmDouble  *recv_buffer  = xNewZeroInit<IssmDouble>(this->Size()*maxdofspernode,"t");
	ISSM_MPI_Request  *send_requests = xNew<ISSM_MPI_Request>(num_procs);
	for (int rank = 0;rank<num_procs;rank++) send_requests[rank] = ISSM_MPI_REQUEST_NULL;

	for(int rank=0;rank<num_procs;rank++){
		if(this->common_send[rank]){
			int  numids = this->common_send[rank];
			send_buffers[rank] = xNew<IssmDouble>(numids*maxdofspernode,"t"); //"t" is required by adolc
			for(int i=0;i<numids;i++){
				int   master_lid = this->common_send_ids[rank][i];
				Node* node=xDynamicCast<Node*>(this->GetObjectByOffset(master_lid));
				_assert_(!node->IsClone());
				for(int j=0;j<node->gsize;j++) send_buffers[rank][i*maxdofspernode+j]=local_ug[node->gdoflist_local[j]];
			}
			ISSM_MPI_Isend(send_buffers[rank],numids*maxdofspernode,ISSM_MPI_DOUBLE,rank,0,IssmComm::GetComm(),&send_requests[rank]);
		}
	}
	for(int rank=0;rank<num_procs;rank++){
		if(this->common_recv[rank]){
			int  numids = this->common_recv[rank];
			ISSM_MPI_Recv(recv_buffer,numids*maxdofspernode,ISSM_MPI_DOUBLE,rank,0,IssmComm::GetComm(),&status);
			for(int i=0;i<numids;i++){
				int   master_lid = this->common_recv_ids[rank][i];
				Node* node=xDynamicCast<Node*>(this->GetObjectByOffset(master_lid));
				for(int j=0;j<node->gsize;j++) local_ug[node->gdoflist_local[j]] = recv_buffer[i*maxdofspernode+j];
			}
		}
	}

	xDelete<IssmDouble>(recv_buffer);
	for(int rank=0;rank<num_procs;rank++){
		if(this->common_send[rank]) ISSM_MPI_Wait(&send_requests[rank],&status);
		xDelete<IssmDouble>(send_buffers[rank]);
	}
	xDelete<IssmDouble*>(send_buffers);
	xDelete<ISSM_MPI_Request>(send_requests);

	/*Assign output pointer*/
	*plocal_ug = local_ug;
}/*}}}*/
