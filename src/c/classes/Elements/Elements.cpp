/*
 * \file Elements.cpp
 * \brief: Implementation of Elements class, derived from DataSet class
 */

/*Headers: {{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./Element.h"
#include "./Elements.h"
#include "../Params/Parameters.h"
#include "../ExternalResults/Results.h"
#include "../ExternalResults/GenericExternalResult.h"
#include "../../toolkits/toolkits.h"
#include "../../shared/shared.h"

using namespace std;
/*}}}*/

/*Object constructors and destructor*/
Elements::Elements(){/*{{{*/
	enum_type=MeshElementsEnum;
	return;
}
/*}}}*/
Elements::~Elements(){/*{{{*/
	return;
}
/*}}}*/

/*Object management*/
void Elements::Configure(Elements* elements,Loads* loads, Nodes* nodes, Vertices* vertices, Materials* materials,Parameters* parameters,Inputs* inputs){/*{{{*/

	vector<Object*>::iterator object;

	for(object=objects.begin() ; object < objects.end(); object++ ){
		Element* element=xDynamicCast<Element*>((*object));
		element->Configure(elements,loads,nodes,vertices,materials,parameters,inputs);
	}

}/*}}}*/
int  Elements::MaxNumNodes(void){/*{{{*/

	int max=0;
	int allmax;
	int numnodes=0;

	/*Now go through all elements, and get how many nodes they own, unless they are clone nodes: */
	for(Object* & object : this->objects){
		Element* element = xDynamicCast<Element*>(object);
		numnodes=element->GetNumberOfNodes();
		if(numnodes>max)max=numnodes;
	}

	/*Grab max of all cpus: */
	ISSM_MPI_Allreduce((void*)&max,(void*)&allmax,1,ISSM_MPI_INT,ISSM_MPI_MAX,IssmComm::GetComm());
	max=allmax;

	return max;
}
/*}}}*/
int  Elements::NumberOfElements(void){/*{{{*/

	int local_nelem;
	int numberofelements;

	local_nelem=this->Size();
	ISSM_MPI_Allreduce((void*)&local_nelem,(void*)&numberofelements,1,ISSM_MPI_INT,ISSM_MPI_SUM,IssmComm::GetComm());

	return numberofelements;
}
/*}}}*/
void Elements::ResetHooks(){/*{{{*/

	vector<Object*>::iterator object;
	Element* element=NULL;

	for ( object=objects.begin() ; object < objects.end(); object++ ){

		element=xDynamicCast<Element*>((*object));
		element->ResetHooks();

	}

}
/*}}}*/
void Elements::SetCurrentConfiguration(Elements* elements,Loads* loads, Nodes* nodes, Vertices* vertices, Materials* materials,Parameters* parameters){/*{{{*/

	vector<Object*>::iterator object;
	Element* element=NULL;

	for ( object=objects.begin() ; object < objects.end(); object++ ){

		element=xDynamicCast<Element*>((*object));
		element->SetCurrentConfiguration(elements,loads,nodes,materials,parameters);

	}

}
/*}}}*/
