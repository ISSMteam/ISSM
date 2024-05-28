/*!\file:  FacesPartitioning.cpp
 * \brief: partition elements and nodes and vertices
 */ 

#include <string.h>
#include "../../classes/classes.h"
#include "../../shared/shared.h"
#include "./ModelProcessorx.h"

void FacesPartitioning(IoModel* iomodel){

	/*If faces are already present, exit*/
	if(iomodel->my_faces) return;

	/*Get faces and elements*/
	CreateFaces(iomodel);
	_assert_(iomodel->elementtofaceconnectivity);

	/*Mesh dependent variables*/
	int elementnbf;
	if(iomodel->domaintype==Domain2DhorizontalEnum){
		elementnbf = 3;
	}
	else if(iomodel->domaintype==Domain2DverticalEnum){
		elementnbf = 3;
	}
	else if(iomodel->domaintype==Domain3DEnum){
		elementnbf = 5;
	}
	else{
		_error_("mesh dimension not supported yet");
	}
	/*output: */
	iomodel->my_faces=xNewZeroInit<bool>(iomodel->numberoffaces);

	for(int i=0;i<iomodel->numberofelements;i++){
		if(iomodel->my_elements[i]){
			for(int j=0;j<elementnbf;j++){
				_assert_(iomodel->elementtofaceconnectivity[i*elementnbf+j] >= 0);
				_assert_(iomodel->elementtofaceconnectivity[i*elementnbf+j] <  iomodel->numberoffaces);
				iomodel->my_faces[iomodel->elementtofaceconnectivity[i*elementnbf+j]] = true;
			}
		}
	}

	if(iomodel->meshelementtype==PentaEnum){
		iomodel->my_vfaces = xNewZeroInit<bool>(iomodel->numberofverticalfaces);
		for(int i=0;i<iomodel->numberofelements;i++){
			if(iomodel->my_elements[i]){
				for(int j=0;j<3;j++){
					_assert_(iomodel->elementtoverticalfaceconnectivity[i*3+j]<iomodel->numberofverticalfaces);
					iomodel->my_vfaces[iomodel->elementtoverticalfaceconnectivity[i*3+j]] = true;
				}
			}
		}
	}
}
