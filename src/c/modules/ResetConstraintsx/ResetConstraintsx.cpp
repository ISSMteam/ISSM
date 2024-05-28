/*!\file ResetConstraintsx
 * \brief: reset thermal penalties
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./ResetConstraintsx.h"
#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"
#include "../ConstraintsStatex/ConstraintsStatex.h"

void ResetConstraintsx(FemModel* femmodel){

	/*Display message*/
	if(VerboseModule()) _printf0_("   Resetting penalties\n");

	/*recover parameters: */
	int analysis_type;
	femmodel->parameters->FindParam(&analysis_type,AnalysisTypeEnum);

	/*Deal with rift first*/
	if(femmodel->loads->numrifts){
		_error_("rift constraints reset not supported yet!");
	}

	/*Reset pengrid to inactive mode*/
	for(Object* & object : femmodel->loads->objects){
		Load* load=(Load*)object;
		if(load->ObjectEnum()==PengridEnum){
			Pengrid* pengrid=(Pengrid*)load;
			pengrid->ResetConstraint();
		}
	}
}
void ResetZigzagCounterx(FemModel* femmodel){

	/*Display message*/
	if(VerboseModule()) _printf0_("   Resetting penalties\n");

	/*Deal with rift first*/
	if(femmodel->loads->numrifts){
		_error_("rift constraints reset not supported yet!");
	}

	/*Reset pengrid to inactive mode*/
	for(Object* & object: femmodel->loads->objects){
		Load* load=(Load*)object;
		if(load->ObjectEnum()==PengridEnum){
			Pengrid* pengrid=(Pengrid*)load;
			pengrid->ResetZigzagCounter();
		}
	}
}
