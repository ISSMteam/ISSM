/*!\file:  CreateSingleNodeToElementConnectivity.cpp
 * \brief: create connectivity table
 */ 

#include "../../shared/shared.h"
#include "../../shared/io/io.h"
#include "../../classes/classes.h"
#include "./ModelProcessorx.h"

void CreateSingleNodeToElementConnectivity(IoModel* iomodel){

	/*Intermediary*/
	int vertexid;
	int elementswidth;

	/*output*/
	int* connectivity=NULL;

	/*Return if connectivity already present*/
	if(iomodel->singlenodetoelementconnectivity) return;

	/*Some checks if debugging*/
	_assert_(iomodel->numberofvertices);
	_assert_(iomodel->numberofelements);
	_assert_(iomodel->my_elements);
	_assert_(iomodel->elements);

	/*Allocate ouput*/
	connectivity=xNewZeroInit<int>(iomodel->numberofvertices);

	/*Get element width*/
	switch(iomodel->meshelementtype){
		case TriaEnum:  elementswidth=3; break;
		case PentaEnum: elementswidth=6; break;
		case TetraEnum: elementswidth=4; break;
		default:  _error_("mesh type "<<EnumToStringx(iomodel->domaintype)<<" not supported yet");
	}

	/*Create connectivity table*/
	for(int i=0;i<iomodel->numberofelements;i++){
		/*!! in parallel we do not want the vertex to be connected to an element that is not in its partition!!*/
		if(iomodel->my_elements[i]){
			for(int j=0;j<elementswidth;j++){
				vertexid=iomodel->elements[elementswidth*i+j];
				_assert_(vertexid>0 && vertexid-1<iomodel->numberofvertices);
				connectivity[vertexid-1]=i+1;
			}
		}
	}

	/*Assign to iomodel*/
	iomodel->singlenodetoelementconnectivity=connectivity;
}
