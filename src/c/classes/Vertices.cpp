/*
 * \file Vertices.cpp
 * \brief: Implementation of Vertices class, derived from DataSet class.
 */

/*Headers: {{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include <vector>
#include <functional>
#include <algorithm>
#include <iostream>

#include "./Vertices.h"
#include "../shared/shared.h"
#include "./Vertex.h"

using namespace std;
/*}}}*/

/*Object constructors and destructor*/
Vertices::Vertices(){/*{{{*/
	this->enum_type              = VerticesEnum;
	this->common_recv            = NULL;
	this->common_recv_ids        = NULL;
	this->common_send            = NULL;
	this->common_send_ids        = NULL;
	this->numberofvertices       = -1;
	this->numberofvertices_local = -1;
	this->numberofmasters_local  = -1;
	return;
}
/*}}}*/
Vertices::~Vertices(){/*{{{*/

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

/*Numerics management*/
Vertices* Vertices::Copy() {/*{{{*/

	int num_proc = IssmComm::GetSize();

	/*Copy dataset*/
	Vertices* output=new Vertices();
	output->sorted    = this->sorted;
	output->numsorted = this->numsorted;
	output->presorted = this->presorted;
	for(vector<Object*>::iterator obj=this->objects.begin() ; obj < this->objects.end(); obj++ ) output->AddObject((*obj)->copy());

	/*Build id_offsets and sorted_ids*/
	output->id_offsets=NULL;
	output->sorted_ids=NULL;
	int objsize = this->numsorted;
	if(this->sorted && objsize>0 && this->id_offsets){
		output->id_offsets=xNew<int>(objsize);
		xMemCpy<int>(output->id_offsets,this->id_offsets,objsize);
	}
	if(this->sorted && objsize>0 && this->sorted_ids){
		output->sorted_ids=xNew<int>(objsize);
		xMemCpy<int>(output->sorted_ids,this->sorted_ids,objsize);
	}

	/*Copy other fields*/
	output->numberofvertices       = this->numberofvertices;
	output->numberofvertices_local = this->numberofvertices_local;
	output->numberofmasters_local  = this->numberofmasters_local;

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
void Vertices::Marshall(MarshallHandle* marshallhandle){ /*{{{*/

	int object_enum = VerticesEnum;
	marshallhandle->call(object_enum);

	/*Check compatibility*/
	int num_procs=IssmComm::GetSize();
	int test = num_procs;
	marshallhandle->call(test);
	if(test!=num_procs) _error_("number of cores is not the same as before");

	marshallhandle->call(this->numberofvertices);
	marshallhandle->call(this->numberofvertices_local);
	marshallhandle->call(this->numberofmasters_local);

	marshallhandle->call(this->common_recv,num_procs);
	marshallhandle->call(this->common_send,num_procs);
	if(marshallhandle->OperationNumber() == MARSHALLING_LOAD){
		this->common_recv_ids = xNew<int*>(num_procs);
		this->common_send_ids = xNew<int*>(num_procs);
	}
	for(int i=0;i<num_procs;i++){
		marshallhandle->call(this->common_recv_ids[i],this->common_recv[i]);
		marshallhandle->call(this->common_send_ids[i],this->common_send[i]);
	}
	DataSet::Marshall(marshallhandle);
}
/*}}}*/
void Vertices::LatLonList(IssmDouble** plat,IssmDouble** plon){/*{{{*/

	/*output: */
	IssmDouble* xyz_serial=NULL;

	/*recover my_rank:*/
	int my_rank=IssmComm::GetRank();

	/*First, figure out number of vertices: */
	int num_vertices=this->NumberOfVertices();

	/*Now, allocate vectors*/
	Vector<IssmDouble>* lat = new Vector<IssmDouble>(num_vertices);
	Vector<IssmDouble>* lon = new Vector<IssmDouble>(num_vertices);

	/*Go through vertices, and for each vertex, object, report it cpu: */
	for(Object* & object : this->objects){
      Vertex* vertex = xDynamicCast<Vertex*>(object);
		lat->SetValue(vertex->sid,vertex->GetLatitude() ,INS_VAL);
		lon->SetValue(vertex->sid,vertex->GetLongitude(),INS_VAL);
	}

	/*Assemble:*/
	lat->Assemble();
	lon->Assemble();

	/*gather on cpu 0: */
	IssmDouble* lat_serial=lat->ToMPISerial();
	IssmDouble* lon_serial=lon->ToMPISerial();

	/*Free resources: */
	*plat = lat_serial;
	*plon = lon_serial;
	delete lat;
	delete lon;
}
/*}}}*/
void Vertices::XYList(IssmDouble** pxcoords,IssmDouble** pycoords){/*{{{*/

	/*output: */
	IssmDouble* xyz_serial=NULL;

	/*recover my_rank:*/
	int my_rank=IssmComm::GetRank();

	/*First, figure out number of vertices: */
	int num_vertices=this->NumberOfVertices();

	/*Now, allocate vectors*/
	Vector<IssmDouble>* xlist = new Vector<IssmDouble>(num_vertices);
	Vector<IssmDouble>* ylist = new Vector<IssmDouble>(num_vertices);

	/*Go through vertices, and for each vertex, object, report it cpu: */
	for(Object* & object : this->objects){
      Vertex* vertex = xDynamicCast<Vertex*>(object);
		xlist->SetValue(vertex->sid,vertex->GetX() ,INS_VAL);
		ylist->SetValue(vertex->sid,vertex->GetY(),INS_VAL);
	}

	/*Assemble:*/
	xlist->Assemble();
	ylist->Assemble();

	/*gather on cpu 0: */
	IssmDouble* x_serial=xlist->ToMPISerial();
	IssmDouble* y_serial=ylist->ToMPISerial();

	/*Free resources: */
	*pxcoords = x_serial;
	*pycoords = y_serial;
	delete xlist;
	delete ylist;
}
/*}}}*/

void Vertices::Finalize(IoModel* iomodel){/*{{{*/

	/*Here we do 3 things:
	 * - count all vertices once for all so that we do not need to call MPI
	 *   every time we need to know the total number of vertices
	 * - Distribute lids (local ids): masters first, slaves second
	 * - Distribute pids (parallel ids)
	 *   */

	/*recover my_rank:*/
	ISSM_MPI_Status status;
	int my_rank   = IssmComm::GetRank();
	int num_procs = IssmComm::GetSize();

	/*1. set number of vertices once for all*/
	this->numberofvertices_local=this->Size();
	this->numberofmasters_local=0;
	for(Object* & object : this->objects){
      Vertex* vertex = xDynamicCast<Vertex*>(object);
		if(!vertex->clone) this->numberofmasters_local++;
	}
	ISSM_MPI_Allreduce((void*)&this->numberofmasters_local,(void*)&this->numberofvertices,1,ISSM_MPI_INT,ISSM_MPI_SUM,IssmComm::GetComm());

	/*2. Distribute lids (First: masters, then clones)*/
	iomodel->my_vertices_lids=xNew<int>(this->numberofvertices);
	for(int i=0;i<this->numberofvertices;i++) iomodel->my_vertices_lids[i] = -1;

	int lid = 0;
	for(Object* & object : this->objects){
      Vertex* vertex = xDynamicCast<Vertex*>(object);
		if(!vertex->clone){
			vertex->lid=lid;
			iomodel->my_vertices_lids[vertex->sid] = lid;
			lid++;
		}
	}
	for(Object* & object : this->objects){
      Vertex* vertex = xDynamicCast<Vertex*>(object);
		if(vertex->clone){
			vertex->lid=lid;
			iomodel->my_vertices_lids[vertex->sid] = lid;
			lid++;
		}
	}

	/*3. Distribute pids based on lids and offsets*/
	int* all_num_masters=xNew<int>(num_procs);
	ISSM_MPI_Gather(&this->numberofmasters_local,1,ISSM_MPI_INT,all_num_masters,1,ISSM_MPI_INT,0,IssmComm::GetComm());
	ISSM_MPI_Bcast(all_num_masters,num_procs,ISSM_MPI_INT,0,IssmComm::GetComm());
	int offset=0;
	for(int i=0;i<my_rank;i++) offset+=all_num_masters[i];
	xDelete<int>(all_num_masters);

	for(Object* & object : this->objects){
      Vertex* vertex = xDynamicCast<Vertex*>(object);
		vertex->pid = vertex->lid+offset;
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
				Vertex* vertex=xDynamicCast<Vertex*>(this->GetObjectByOffset(this->common_send_ids[rank][i]));
				send_truepids[rank][i] = vertex->pid;
			}
			ISSM_MPI_Isend(send_truepids[rank],numids,ISSM_MPI_INT,rank,0,IssmComm::GetComm(),&send_requests[rank]);
		}
	}
	for(int rank=0;rank<num_procs;rank++){
		if(this->common_recv[rank]){
			int  numids = this->common_recv[rank];
			ISSM_MPI_Recv(recv_truepids,numids,ISSM_MPI_INT,rank,0,IssmComm::GetComm(),&status);
			for(int i=0;i<numids;i++){
				Vertex* vertex=xDynamicCast<Vertex*>(this->GetObjectByOffset(this->common_recv_ids[rank][i]));
				vertex->pid = recv_truepids[i];
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
}/*}}}*/
int Vertices::NumberOfVertices(){/*{{{*/
	return this->numberofvertices;
}/*}}}*/
int Vertices::NumberOfVerticesLocal(void){/*{{{*/
	return this->numberofmasters_local;
}/*}}}*/
int Vertices::NumberOfVerticesLocalAll(void){/*{{{*/
	return this->numberofvertices_local;
}/*}}}*/
