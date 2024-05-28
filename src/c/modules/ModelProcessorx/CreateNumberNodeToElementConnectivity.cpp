/*!\file:  CreateNumberNodeToElementConnectivity.cpp
 * \brief: create connectivity table
 */ 

#ifdef HAVE_CONFIG_H
#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../../shared/shared.h"
#include "../../classes/classes.h"
#include "../../shared/io/io.h"
#include "./ModelProcessorx.h"

void CreateNumberNodeToElementConnectivity(IoModel* iomodel){

	/*Intermediary*/
	int i,j;
	int vertexid;
	int elementswidth;

	/*output*/
	int* connectivity=NULL;

	/*Check that this has not been done yet*/
	if(iomodel->numbernodetoelementconnectivity) return;

	/*Some checks if debugging*/
	_assert_(iomodel->numberofvertices);
	_assert_(iomodel->numberofelements);
	_assert_(iomodel->elements);

	/*Allocate ouput*/
	connectivity=xNewZeroInit<int>(iomodel->numberofvertices);

	/*Get element width*/
	switch(iomodel->meshelementtype){
		case TriaEnum:  elementswidth=3; break;
		case TetraEnum: elementswidth=4; break;
		case PentaEnum: elementswidth=6; break;
		default:                   _error_("mesh not supported yet");
	}

	/*Create connectivity table*/
	for (i=0;i<iomodel->numberofelements;i++){
		for (j=0;j<elementswidth;j++){
			vertexid=iomodel->elements[elementswidth*i+j];
			_assert_(vertexid>0 && vertexid-1<iomodel->numberofvertices);
			connectivity[vertexid-1]+=1;
		}
	}

	/*Assign to iomodel*/
	iomodel->numbernodetoelementconnectivity=connectivity;
}
