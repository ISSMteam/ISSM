/*!\file ResetFSBasalBoundaryConditionx
 * \brief: reset coordinate system for full-FS: tangential to the bedrock
 */

#include "./ResetFSBasalBoundaryConditionx.h"
#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"

void ResetFSBasalBoundaryConditionx(Elements* elements,Nodes* nodes, Vertices* vertices, Loads* loads,Materials* materials,Parameters* parameters){

	Element *element = NULL;

for(Object* & object : elements->objects){
      element = xDynamicCast<Element*>(object);
		element->ResetFSBasalBoundaryCondition();
	}

}
