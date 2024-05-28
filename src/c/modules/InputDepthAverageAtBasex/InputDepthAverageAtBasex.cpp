/*!\file InputDepthAverageAtBasex
 * \brief: extrude input
 */

#include "./InputDepthAverageAtBasex.h"
#include "../../shared/shared.h"
#include "../../classes/classes.h"
#include "../../toolkits/toolkits.h"

void InputDepthAverageAtBasex(FemModel* femmodel,int original_enum, int new_enum){
	for(Object* & object : femmodel->elements->objects){
		Element* element=xDynamicCast<Element*>(object);
		element->InputDepthAverageAtBase(original_enum,new_enum);
	}
}
