/*!\file UpdateDynamicConstraintsx
 * \brief module to update single point constraints  out of new spc vector, for next time step.
 */

#include "./UpdateDynamicConstraintsx.h"

#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"

void UpdateDynamicConstraintsx(Constraints* constraints,Nodes* nodes,Parameters* parameters,Vector<IssmDouble>* yg){

	int configuration_type;
	IssmDouble* yg_serial=NULL;

	/*Get current configuration*/
	parameters->FindParam(&configuration_type,ConfigurationTypeEnum);

	/*serialize yg, so nodes can index into it: */
	yg_serial=yg->ToMPISerial();

	for(Object* & object : constraints->objects){
		Constraint* constraint=(Constraint*)object;
		if(constraint->ObjectEnum()==SpcDynamicEnum){
			((SpcDynamic*)constraint)->SetDynamicConstraint(nodes,yg_serial);
		}

	}

	/*Free resources:*/
	xDelete<IssmDouble>(yg_serial);
}
