/*!\file InputUpdateFromMatrixDakotax
 * \brief: update datasets using  parameter inputs
 */

#include "./InputUpdateFromMatrixDakotax.h"
#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"
#include "../InputUpdateFromVectorDakotax/InputUpdateFromVectorDakotax.h"

void InputUpdateFromMatrixDakotax(FemModel* femmodel,double* matrix,int nrows,int ncols, int name, int type){

	int numberofvertices,numberofelements;

	numberofvertices=femmodel->vertices->NumberOfVertices();
	numberofelements=femmodel->elements->NumberOfElements();

	if((ncols==1) && (nrows==numberofvertices || nrows==numberofelements)) InputUpdateFromVectorDakotax(femmodel,matrix,name,type);
	else{

		/*Update elements, nodes, loads and materials from inputs: */
		for(Object* & object : femmodel->elements->objects){
			Element* element=xDynamicCast<Element*>(object);
			element->InputUpdateFromMatrixDakota(matrix,nrows,ncols,name,type);
		}
	}
}
