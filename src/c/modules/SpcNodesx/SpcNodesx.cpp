/*!\file SpcNodesx
 * \brief: establish single point constraints on all nodes, as well as constraints vector.
 */

#include "./SpcNodesx.h"

#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"

void SpcNodesx(Nodes* nodes,Constraints* constraints,Parameters* parameters){

	for(Object* & object: constraints->objects){
		Constraint* constraint=xDynamicCast<Constraint*>(object);
		constraint->ConstrainNode(nodes,parameters);
	}
}
