/*
 * \file Loads.cpp
 * \brief: Implementation of Loads class, derived from DataSet class.
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

#include "../../shared/io/Comm/IssmComm.h"
#include "../../shared/Numerics/recast.h"
#include "../../shared/Enum/EnumDefinitions.h"
#include "../../shared/Exceptions/exceptions.h"
#include "../../shared/MemOps/MemOps.h"
#include "../../shared/io/Marshalling/Marshalling.h"
#include "./Loads.h"
#include "./Load.h"

using namespace std;
/*}}}*/

/*Object constructors and destructor*/
Loads::Loads(){/*{{{*/
	this->enum_type=LoadsEnum;
	this->numrifts     = 0;
	this->numpenalties = 0;
	return;
}
/*}}}*/
Loads::~Loads(){/*{{{*/
	return;
}
/*}}}*/

Loads* Loads::Copy() {/*{{{*/

	int num_proc = IssmComm::GetSize();

	/*Copy dataset*/
	Loads* output=new Loads();
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
	output->numrifts = this->numrifts;
	output->numpenalties = this->numpenalties;

	return output;
}
/*}}}*/
void  Loads::Marshall(MarshallHandle* marshallhandle){ /*{{{*/

	int object_enum = LoadsEnum;
	marshallhandle->call(object_enum);
	marshallhandle->call(this->numrifts);
	marshallhandle->call(this->numpenalties);

	DataSet::Marshall(marshallhandle);
}
/*}}}*/

/*Numerics:*/
void Loads::Configure(Elements* elements,Loads* loads, Nodes* nodes, Vertices* vertices, Materials* materials,Parameters* parameters){/*{{{*/

	vector<Object*>::iterator object;
	for(object=objects.begin() ; object < objects.end(); object++){
		Load* load=xDynamicCast<Load*>(*object);
		load->Configure(elements,loads,nodes,vertices,materials,parameters);
	}
}
/*}}}*/
void Loads::Finalize(){/*{{{*/

	/*Count Rifts and penalties*/
	int ispenalty=0;
	int isrift=0;
	int allcount;

	/*Now go through all loads, and get how many nodes they own, unless they are clone nodes: */
	for(Object* & object : this->objects){
      Load* load = xDynamicCast<Load*>(object);
		if(load->IsPenalty()){
			ispenalty++;
		}
      if(load->ObjectEnum()==RiftfrontEnum){
         isrift++;
      }
	}

	/*Grab sum of all cpus: */
	ISSM_MPI_Allreduce((void*)&ispenalty,(void*)&allcount,1,ISSM_MPI_INT,ISSM_MPI_SUM,IssmComm::GetComm());
	this->numpenalties = allcount;

	ISSM_MPI_Allreduce((void*)&isrift,(void*)&allcount,1,ISSM_MPI_INT,ISSM_MPI_SUM,IssmComm::GetComm());
	this->numrifts= allcount;

}
/*}}}*/
bool Loads::IsPenalty(){/*{{{*/

	if(this->numpenalties>0)
	 return true;
	else
	 return false;
}
/*}}}*/
int  Loads::MaxNumNodes(){/*{{{*/

	int max=0;
	int allmax;

	/*Now go through all loads, and get how many nodes they own, unless they are clone nodes: */
	for(Object* & object : this->objects){
      Load* load = xDynamicCast<Load*>(object);
		int numnodes=load->GetNumberOfNodes();
		if(numnodes>max)max=numnodes;
	}

	/*Grab max of all cpus: */
	ISSM_MPI_Allreduce((void*)&max,(void*)&allmax,1,ISSM_MPI_INT,ISSM_MPI_MAX,IssmComm::GetComm());
	return allmax;
}
/*}}}*/
int  Loads::NumberOfLoads(void){/*{{{*/

	int localloads;
	int numberofloads;

	/*Get number of local loads*/
	localloads=this->Size();

	/*figure out total number of loads combining all the cpus (no clones here)*/
	ISSM_MPI_Reduce(&localloads,&numberofloads,1,ISSM_MPI_INT,ISSM_MPI_SUM,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&numberofloads,1,ISSM_MPI_INT,0,IssmComm::GetComm());

	return numberofloads;
}
/*}}}*/
void Loads::ResetHooks(){/*{{{*/

	vector<Object*>::iterator object;
	Load* load=NULL;

	for ( object=objects.begin() ; object < objects.end(); object++ ){

		load=xDynamicCast<Load*>((*object));
		load->ResetHooks();

	}

}
/*}}}*/
void Loads::SetCurrentConfiguration(Elements* elements,Loads* loads, Nodes* nodes, Vertices* vertices, Materials* materials,Parameters* parameters){/*{{{*/

	vector<Object*>::iterator object;
	Load* load=NULL;

	for ( object=objects.begin() ; object < objects.end(); object++ ){
		load=xDynamicCast<Load*>(*object);
		load->SetCurrentConfiguration(elements,loads,nodes,vertices,materials,parameters);
	}
}
/*}}}*/
