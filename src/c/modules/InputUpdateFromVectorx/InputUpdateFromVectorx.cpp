/*!\file InputUpdateFromVectorx
 * \brief: update datasets using  parameter inputs
 */

#include "./InputUpdateFromVectorx.h"
#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"

void InputUpdateFromVectorx(FemModel* femmodel,Vector<IssmDouble>* vector, int name, int type){

	if(type==VertexPIdEnum){
		IssmDouble* serial_vector=NULL;
		femmodel->GetLocalVectorWithClonesVertices(&serial_vector,vector);
		InputUpdateFromVectorx(femmodel,serial_vector,name,VertexLIdEnum);
		xDelete<IssmDouble>(serial_vector);
	}
	else{
		IssmDouble* serial_vector=vector->ToMPISerial();
		InputUpdateFromVectorx(femmodel,serial_vector,name,type);
		xDelete<IssmDouble>(serial_vector);
	}
}

void InputUpdateFromVectorx(FemModel* femmodel,IssmDouble* vector, int name, int type){

	/*Update elements, nodes, loads and materials from inputs: */
	for(Object* & object : femmodel->elements->objects){
		Element* element=xDynamicCast<Element*>(object);
		element->InputUpdateFromVector(vector,name,type);
	}
}
