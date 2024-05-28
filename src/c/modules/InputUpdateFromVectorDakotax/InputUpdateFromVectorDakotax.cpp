/*!\file InputUpdateFromVectorDakotax
 * \brief: update datasets using  parameter inputs
 */

#include "./InputUpdateFromVectorDakotax.h"
#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"

void InputUpdateFromVectorDakotax(FemModel* femmodel,Vector<IssmDouble>* vector, int name, int type){

	IssmDouble* serial_vector=vector->ToMPISerial();
	InputUpdateFromVectorDakotax(femmodel,serial_vector,name, type);

	/*Free resources:*/
	xDelete<double>(serial_vector);
}

void InputUpdateFromVectorDakotax(FemModel* femmodel,IssmDouble* vector, int name, int type){

	/*Update elements, nodes, loads and materials from inputs: */
	for(Object* & object : femmodel->elements->objects){
		Element* element=xDynamicCast<Element*>(object);
		element->InputUpdateFromVectorDakota(vector,name,type);
	}
}
