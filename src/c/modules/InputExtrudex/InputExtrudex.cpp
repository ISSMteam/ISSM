/*!\file InputExtrudex
 * \brief: extrude input
 */

#include "./InputExtrudex.h"
#include "../../shared/shared.h"
#include "../../classes/classes.h"
#include "../../toolkits/toolkits.h"

void InputExtrudex(FemModel* femmodel,int input_enum,int start){
	for(Object* & object : femmodel->elements->objects){
		Element* element=xDynamicCast<Element*>(object);
		element->InputExtrude(input_enum,start);
	}
}
