/*!\file:  EdgesPartitioning.cpp
 * \brief: partition elements and nodes and vertices
 */ 

#include <string.h>
#include "../../classes/classes.h"
#include "../../shared/shared.h"
#include "./ModelProcessorx.h"

void EdgesPartitioning(IoModel* iomodel){

	/*If faces are already present, exit*/
	if(iomodel->my_edges) return;

	/*Get edges and elements*/
	CreateEdges(iomodel);
	_assert_(iomodel->elementtoedgeconnectivity);

	/*Mesh dependent variables*/
	int elementnbe;
	switch(iomodel->meshelementtype){
		case TriaEnum:  elementnbe = 3; break;
		case TetraEnum: elementnbe = 6; break;
		case PentaEnum: elementnbe = 9; break;
		default: _error_("mesh dimension not supported yet");
	}

	/*output: */
	iomodel->my_edges  = xNewZeroInit<bool>(iomodel->numberofedges);

	for(int i=0;i<iomodel->numberofelements;i++){
		if(iomodel->my_elements[i]){
			for(int j=0;j<elementnbe;j++){
				iomodel->my_edges[iomodel->elementtoedgeconnectivity[i*elementnbe+j]] = true;
			}
		}
	}

	if(iomodel->meshelementtype==PentaEnum){
		iomodel->my_vedges = xNewZeroInit<bool>(iomodel->numberofverticaledges);
		iomodel->my_hedges = xNewZeroInit<bool>(iomodel->numberofhorizontaledges);
		for(int i=0;i<iomodel->numberofelements;i++){
			if(iomodel->my_elements[i]){
				for(int j=0;j<3;j++) iomodel->my_vedges[iomodel->elementtoverticaledgeconnectivity[i*3+j]]   = true;
				for(int j=0;j<6;j++) iomodel->my_hedges[iomodel->elementtohorizontaledgeconnectivity[i*6+j]] = true;
			}
		}
	}
}
